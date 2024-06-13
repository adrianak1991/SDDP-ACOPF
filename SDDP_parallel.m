%%% Here we implement SDDP algorithm for the multistage stochastic optimization
%%% problem including OPF component as defined in the paper: "Stochastic
%%% Dual Dynami Programming For Optimal Power Flow Problems under
%%% Uncertainty"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic; clear; clc; yalmip('clear'); 

% set random seed for reproducability of serial version of code
param.seed = 53452; rng(param.seed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 1. Set Parameters                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
area = 1;                               % area of network to use for the problem
param.D = 7;                            % number of random stages
param.H = 24;                           % number of substages
param.N = 100;                          % number of nodes per (non-initial) stage in lattice
param.num_iter = 100;                   % number of iterations of the SDDP algorithm
param.ub_scen = 150;                    % number of forward passes for calculation of ub
param.ub_freq = 5;                      % frequency of ub calculation
param.EV_4D = false;                    % boolean that indicates whether a valid voltage vector should be searched for in the 4 dimensional null space
param.day_0 = 191;                      % day of year which is starting point for the scenario lattice (one day before that)
param.benchmark = 0.95;                 % benchmark from PCA to choose proper number of variables
param.num_scen = 1e5;                   % number of samples used to estimate scenario lattice
param.curtailment_cost_plus = 1000 * 100; % Curtailment costs (demand)
param.curtailment_cost_minus = 100 * 100; % Curtailment costs (generation)
param.cost_storage = 230 * 100;         % storage cost per MWh capacity
param.B_0 = 0.75;                       % initial storage level
param.eta = 1;                          % charge/discharge efficiency
param.new_fit = false;                  % whether the econometric model should be fitted or loaded - quasar required
param.VRES_mult = 2;                    % how many times increase of renewables

param.log_file_name = sprintf("%d-%02.0f-%02.0f %02.0f-%02.0f-%02.0f.log", year(now), month(now), day(now), hour(now), minute(now), second(now));
param.log_file = fopen(param.log_file_name, 'w');

% Solver parameters
acc = 10^(-11);
param.options = sdpsettings('solver','mosek','verbose',0, 'cachesolvers', 1,'MOSEK.MSK_DPAR_INTPNT_CO_TOL_REL_GAP',acc,'MOSEK.MSK_DPAR_INTPNT_CO_TOL_DFEAS',acc,...
    'MOSEK.MSK_DPAR_INTPNT_CO_TOL_INFEAS',acc,'MOSEK.MSK_DPAR_INTPNT_CO_TOL_MU_RED', acc,'MOSEK.MSK_DPAR_INTPNT_CO_TOL_PFEAS',acc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 2. Get Data, Estimate Model, Learn Lattice %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get data
data = prepareData(area);

% Create lattice from simulations from a model based on historical data
[data_scenarios, probability_transitions,N_final] = createLattice(data, param);
param.N_final = N_final;
[Pd_scen,Qd_scen] = calculateDemand(data, param, data_scenarios);

data.prob = probability_transitions;
data.Pd_scen = Pd_scen;
data.Qd_scen = Qd_scen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 3. Set Up Problems                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = data.net.k; % number of buses with storage

% parameters of cuts for first stage (d=0) where we decide about storage size
e0_0 = [];                      % intercepts
e1_0 = [];                      % slopes wrt storage capacity
sol_cap_0 = [];                 % storage decisions

% parameters of cuts for resource state (d = 1..6)
e0 = cell(param.N, param.D);          % intercepts
e1 = cell(param.N, param.D);          % slopes wrt storage capacity (for all the nodes)
e2 = cell(param.N, param.D);          % slopes wrt to storage levels
sol_level = zeros(param.num_iter, k, param.D + 1);

% bounds variables
forward_val_vec = cell(param.num_iter,1);
upper_bound = zeros(param.num_iter,1);
lower_bound = zeros(param.num_iter,1);
confidence_ub = zeros(param.num_iter,2);

% status variables
forward_status = cell(param.num_iter, 1);
backward_status = cell(param.num_iter, 1);

% Assign scenarios to workers for forward and backward passes
my_cluster = gcp;
param.num_workers = my_cluster.NumWorkers;
assignment_bw = cell(param.num_workers, param.D);
assignment_fw = cell(param.num_workers, 1);

s_bw = floor(param.N_final/param.num_workers); s_fw = floor(param.ub_scen/param.num_workers);
for w = 1:param.num_workers
    for stage=1:param.D
        assignment_bw{w,stage} = s_bw(stage)*(w-1)+1:s_bw(stage)*w;
        % distribute the remaining scenarios to the first mod(param.N, param.num_workers)
        % workers
        if w <= mod(param.N_final(stage), param.num_workers)
            assignment_bw{w,stage} = [assignment_bw{w,stage}, param.num_workers*s_bw(stage)+w];
        end
    end

    assignment_fw{w} = s_fw*(w-1)+1:s_fw*w;
    if w <= mod(param.ub_scen, param.num_workers)
        assignment_fw{w} = [assignment_fw{w}, param.num_workers*s_fw+w];
    end
end

yalmiptime = 0; solvertime = 0;
iterationtime = zeros(param.num_iter,1);

fprintf('%-5s%-8s%-11s%-11s%-10s%-9s%-9s\n', 'i', 'time', 'storage', 'lower', 'upper', 'upper_l', 'upper_u');
fprintf('-------------------------------------------------------------\n');
for iter = 1:param.num_iter
    %%%%% FORWARD PASS
    % First solve the first stage problem deciding about storage capacities
    [storage_size, first_stage_cost, first_stage_objective, time] = firstStage(data, param, e0_0, e1_0);
    yalmiptime = yalmiptime + time.yalmiptime;
    solvertime = solvertime + time.solvertime;
    sol_cap_0(iter,:) = storage_size;
    lower_bound(iter) = first_stage_objective;

    % distinguish between the iterations with ub calculations and those without
    if mod(iter, param.ub_freq) == 0
        % parallel forward pass
        spmd
            [levels, forward_val, status, nodes, storage_hourly,curtailment_hourly] = forwardPass_parallel(data, param, length(assignment_fw{labindex}), sol_cap_0(iter,:)', e0, e1, e2);
        end
        levels_result = cat(3,levels{:}) ;
        sol_level(iter,:,:) = levels_result(:,:,1);
        forward_val_vec{iter} = first_stage_cost + [forward_val{:}];
        status_result = [status{:}];
        forward_status{iter} = [status_result.num_errors];
        yalmiptime = yalmiptime + max([status_result.yalmiptime]);
        solvertime = solvertime + max([status_result.solvertime]);

        nodes_result = [nodes{:}];
        storage_hourly_result =  cat(4,storage_hourly{:}) ;
        curtailment_hourly_result = [curtailment_hourly{:}];

        % Calculation of upper bound based on the forward pass
        upper_bound(iter) = mean(forward_val_vec{iter});
        confidence_ub(iter, :) = upper_bound(iter) - 1.96 * std(forward_val_vec{iter})/sqrt(param.ub_scen).*[1 -1];
    else
        [levels, forward_val, status] = forwardPass_parallel(data, param, 1, sol_cap_0(iter,:)', e0, e1, e2);
        sol_level(iter,:,:) = levels(:,:,1);
        forward_val_vec{iter} = first_stage_cost + forward_val;
        forward_status{iter} = status.num_errors;
        yalmiptime = yalmiptime + status.yalmiptime;
        solvertime = solvertime + status.solvertime;

        % Calculation of upper bound based on the forward pass
        upper_bound(iter) = forward_val_vec{iter};
    end

    %% BACKWARD PASS
    for time=param.D:-1:1
        spmd
            [e0_val,e1_val,e2_val,status] = backwardPass_parallel(data, param, assignment_bw{labindex,time}', time, sol_cap_0(iter,:)', sol_level(iter,:,time,1), e0, e1, e2);
        end
        order_bw = [assignment_bw{:,time}];
        single_e0(1,order_bw) = [e0_val{:}];
        single_e1(:,order_bw)= [e1_val{:}];
        single_e2(:,order_bw) = [e2_val{:}];

        status_result = [status{:}];
        backward_status{iter} = [backward_status{iter}; sum([status_result.num_errors])];
        yalmiptime = yalmiptime + max([status_result.yalmiptime]);
        solvertime = solvertime + max([status_result.solvertime]);

        % Getting a slope and intercept of cut for given stage (expected
        % value over scenarios)
        prob_matrix = probability_transitions{time,1};
        if (time>1)
            for n = 1:param.N_final(time-1)
                e0{n, time-1} = [e0{n, time-1}; prob_matrix(n,:) * single_e0(1,1:param.N_final(time))'];
                e1{n, time-1} = [e1{n, time-1}; prob_matrix(n,:) * single_e1(:,1:param.N_final(time))'];
                e2{n, time-1} = [e2{n, time-1}; prob_matrix(n,:) * single_e2(:,1:param.N_final(time))'];
            end
        else
            e0_0(iter,1) =  (prob_matrix * single_e0(1,1:param.N_final(time))')';
            e1_0(iter,:) =  (prob_matrix * single_e1(:,1:param.N_final(time))')' + param.B_0 * (prob_matrix * single_e2(:,1:param.N_final(time))')';
        end
    end

    if mod(iter, param.ub_freq) == 0
        fprintf('%-3d  %-6.2f  %-9.2f  %-9.0f  %-8.0f  %-7.0f  %-7.0f\n', iter, toc/60, sum(sol_cap_0(iter,:)) * 100, lower_bound(iter), upper_bound(iter), confidence_ub(iter, 1), confidence_ub(iter, 2));
    else
        fprintf('%-3d  %-6.2f  %-9.2f  %-9.0f  %-8.0f\n', iter, toc/60, sum(sol_cap_0(iter,:)) * 100, lower_bound(iter), forward_val_vec{iter});
    end
    iterationtime(iter) = toc;

    fprintf('Fraction of problems in backward pass: %.2f\n', sum(backward_status{iter}, 'all') / (sum(param.N_final) * param.H));

    fprintf(param.log_file, '----------------------------- Finisihed iteration %d. -----------------------------\n', iter);
end


%% Plots
% - plot upper and lower bounds
I = (2*param.ub_freq:param.ub_freq:param.num_iter)';
xCoord = [I; flipud(I)];
inBetween = [confidence_ub(I, 1); flipud(confidence_ub(I, 2))];
figure; plot(I, [upper_bound(I), lower_bound(I)]); hold on;
h = fill(xCoord, inBetween, 'b', 'LineStyle', 'none'); set(h,'facealpha',.3); hold on;
legend('upper bound', 'lower bound');title('Convergence Plot');xlim([I(1),I(end)]);
set(gca,'FontSize',16)

% - plot number of problems with rank conditions in stages for forward and
% backward passes
prob_fw = zeros(param.num_iter, 1);
prob_bw = zeros(param.num_iter, 1);
for i = 1:param.num_iter
    if mod(i, param.ub_freq) == 0
        prob_fw(i) = sum(forward_status{i}, 'all') / (param.H * param.D * param.ub_scen);
    else
        prob_fw(i) = sum(forward_status{i}, 'all') / (param.H * param.D);
    end
    prob_bw(i) = sum(backward_status{i}, 'all') / (sum(param.N_final) * param.H);
end
figure; plot([prob_fw, prob_bw])
legend('Problems Forward Pass', 'Problems Backward Pass');

% - plot storage decisions over time (from the last convergence check) and
% the associated power levels at the node
storage = data.net.k;
storage_level = reshape(permute(squeeze(storage_hourly_result(storage, :, :, :)),[2,1,3]), param.H * param.D, param.ub_scen);
figure; subplot(2,1,1); fanChart(1:param.D*param.H, storage_level*100); title('Storage Levels');xlim([0,168]);
xlabel("Time (h)"); ylabel("Storage level (MWh)");
set(gca,'FontSize',16)

node = data.net.storage_idx(storage);
residual_demand = zeros(param.H * param.D, param.ub_scen);
for scen = 1:param.ub_scen
    for stage = 1:param.D
        residual_demand((stage-1)*24+1:stage*24, scen) = squeeze(data.Pd_scen{stage}(nodes_result(stage, scen), node, :));
    end
end
subplot(2,1,2); fanChart(1:param.D*param.H, residual_demand*100); title('Residual Active Power Demand');xlim([0,168]);
xlabel("Time (h)"); ylabel("Demand (MW)");
set(gca,'FontSize',16)

% - plot curtailment values over time (from the last convergence check)
bus_id = 1; % choose value between 1 and data.net.n
ylim_values = [0,10^-6];
curtailment_p_plus = cat(4,curtailment_hourly_result.p_plus);
curtailment_p_plus_level = reshape(squeeze(curtailment_p_plus(bus_id, :, :, :)), param.H * param.D, param.ub_scen);
figure; subplot(3,1,1); fanChart(1:param.D*param.H, curtailment_p_plus_level); title('Demand Curtailment Levels'); ylim(ylim_values);

curtailment_p_minus = cat(4,curtailment_hourly_result.p_minus);
curtailment_p_minus_level = reshape(squeeze(curtailment_p_minus(bus_id, :, :, :)), param.H * param.D, param.ub_scen);
subplot(3,1,2); fanChart(1:param.D*param.H, curtailment_p_minus_level); title('Generation Curtailment Levels'); ylim(ylim_values);

residual_demand = zeros(param.H * param.D, param.ub_scen);
for scen = 1:param.ub_scen
    for stage = 1:param.D
        residual_demand((stage-1)*24+1:stage*24, scen) = squeeze(data.Pd_scen{stage}(nodes_result(stage, scen), bus_id, :));
    end
end
subplot(3,1,3); fanChart(1:param.D*param.H, residual_demand); title('Residual Active Power Demand');

toc

%% Save data
save(sprintf("run_%d-%02.0f-%02.0f %02.0f-%02.0f.mat", year(now), month(now), day(now), hour(now), minute(now)));
