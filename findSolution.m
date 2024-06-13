%%% This function solves the problem for the given (number of) scenarios
%%% and provide the solution of the investment planning (if not provided) and
%%% for the operational planning.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cost, W, P, Q, P_g, Q_g, B_k, curtailment, B_max, nodes, status] = findSolution(data, param, num_scen, e0_0, e1_0, e0, e1, e2, B_max_value, nodes_value)

% loading all network related data into the local workspace of the function
struct2vars(data.net);
struct2vars(data.admittance);

param.log_file = fopen(param.log_file_name, 'w');

% define base components of optimization problem
[cons_bas,h_bas, dec] = defineProblem(data, param);

% number of cuts presently in the VFA
num_cuts = length(e0{1,1});

% status information
status.errorCodes = ones(param.D,param.H,num_scen); % Here we collect information if the condition for matrix A is satisfied in the given step
status.num_errors = zeros(param.D);
status.solvertime = 0; status.yalmiptime = 0;


W = cell(param.H,param.D,num_scen);
P = zeros(n,param.H,param.D,num_scen);
Q = zeros(n,param.H,param.D,num_scen);
P_g = zeros(m,param.H, param.D,num_scen);
Q_g = zeros(m,param.H, param.D,num_scen);
B_k = zeros(k,param.H,param.D,num_scen);

% hourly curtailment values
curtailment = cell(num_scen,1);
ub_curtailment_cost = 0;

%% First solve the first stage problem deciding about storage capacities
if ~exist('B_max_value','var')
    [B_max,first_stage_cost,~,time] = firstStage(data, param, e0_0, e1_0);
    status.yalmiptime = status.yalmiptime + time.yalmiptime;
    status.solvertime = status.solvertime + time.solvertime;

    cost = first_stage_cost*ones(num_scen,1);
else
    B_max = B_max_value;
    cost = param.cost_storage*sum(B_max)*ones(num_scen,1);
end

% Resource states for all scenarios (initial levels at stage). since we also
% want to store the storage level of at the end of planning, i.e., after stage
% param.D, we need to store the initial level of stage param.D+1.
new_sol = zeros(k, param.D+1,num_scen);
new_sol(:,1,:) = repmat(param.B_0*B_max,1,1,num_scen);

%% GENERATING SCENARIOS IF NOT PROVIDED
if ~exist('nodes_value','var')
    nodes = ones(param.D,num_scen);
    for stage = 1:param.D
        for i=1:num_scen
                % simulate a node for the stage
                if stage == 1
                    nodes(stage,i) = find(cumsum(data.prob{stage}(1,:)) >= rand, 1);
                else
                    nodes(stage,i) = find(cumsum(data.prob{stage}(nodes(stage-1,i),:)) >= rand, 1);
                end
        end
    end
else
    nodes = nodes_value;
end

%% PARAMETERS SPECIFIC FOR THIS PROBLEM
h_0 = h_bas - sum((dec.kappa_u + dec.beta_u_minus * param.eta/4 + dec.beta_u_plus/(4*param.eta))' * B_max);

if num_cuts > 0
    % Lagrange multiplier associated with value function approximation
    % constraint.
    nu = sdpvar(num_cuts,1);
    cons_0 = [cons_bas, nu>=0, (sum(nu)==1):'VFA'];
else
    cons_0 = cons_bas;
end

%% RUNNING OPTIMIZATION PROBLEM FOR DIFFERENT SCENARIOS
for stage = 1:param.D
    for i=1:num_scen
        % We add part of objective function involving the random parameter
        h = h_0 - dec.sigma(:,1)' * new_sol(:, stage,i);
        h = h + sum(dec.lambda .* squeeze(data.Pd_scen{stage}(nodes(stage,i),:,:)),'all')+sum(dec.gamma .* squeeze(data.Qd_scen{stage}(nodes(stage,i),:,:)),'all');

        % if there is already a VFA and its not the last stage (where there
        % is no VFA per definition)
        if num_cuts > 0 && stage < param.D
            h = h + e0{nodes(stage,i),stage}' * nu;
            h = h + nu' * e1{nodes(stage,i),stage} * B_max;

            % constraint with dual multiplier of last storage level
            cons = [cons_0, (e2{nodes(stage,i), stage}' * nu + dec.sigma(:, param.H) + dec.kappa_u(:, param.H) >= 0):sprintf('storage_%d', param.H)];
        else
            if stage==param.D
                tau = sdpvar(k,1,'full');
    
                h = h - tau'*(param.B_0*B_max);
                cons = [cons_0, (dec.sigma(:, param.H) + dec.kappa_u(:, param.H) + tau >= 0):sprintf('storage_%d', param.H)];
            else
                % constraint with dual multiplier of last stage's storage level
                % if there is no VFA or it is the last stage
                cons = [cons_0, (dec.sigma(:, param.H) + dec.kappa_u(:, param.H) >= 0):sprintf('storage_%d', param.H)];
            end
        end

        diagnostics = optimize(cons, -h, param.options);
        status.solvertime = status.solvertime + diagnostics.solvertime;
        status.yalmiptime = status.yalmiptime + diagnostics.yalmiptime;

        new_acc = -inf; new_options = param.options;
        while (diagnostics.problem ~= 0) && new_acc < 10^(-9)
            % decrease numerical accuracy
            new_acc = new_options.mosek.MSK_DPAR_INTPNT_CO_TOL_INFEAS * 10;
            new_options = sdpsettings('solver','mosek','verbose',0,'MOSEK.MSK_DPAR_INTPNT_CO_TOL_REL_GAP',new_acc,'MOSEK.MSK_DPAR_INTPNT_CO_TOL_DFEAS',new_acc,...
                'MOSEK.MSK_DPAR_INTPNT_CO_TOL_INFEAS',new_acc,'MOSEK.MSK_DPAR_INTPNT_CO_TOL_MU_RED', new_acc,'MOSEK.MSK_DPAR_INTPNT_CO_TOL_PFEAS',new_acc);
            diagnostics = optimize(cons, -h, new_options);
        end

        if(diagnostics.problem ~= 0)
            error("F: Solver error: " + diagnostics.info);
        else
            % We have to check if the condition on the matrix A
            % is satisfied (null space has dimension 2).
            for hour = 1:param.H
                A_sol = value(dec.matA(:,:,hour));
                [eig_vec, eig_val] = eig(A_sol);

                W = dual(cons(sprintf('W_%d', hour)));
                [count, max_u, max_l, mean_dev, max_P, max_Q, V, curtailment_cost] = checkPrimal_power(data, param, eig_vec, double(dec.mu_u(:,hour)), double(dec.mu_l(:,hour)), W);
                ub_curtailment_cost = ub_curtailment_cost + curtailment_cost;

                if max_u > 1e-3 || max_l > 1e-3 || max_P > 1e-3 || max_Q > 1e-3
                  status.errorCodes(stage,hour,i) = 1;
                  status.num_errors(stage) = status.num_errors(stage) + 1;
                end

                % Transformation into optimal powers
                P(:,hour,stage,i) = real(V .* (conj(Y) * conj(V)));
                Q(:,hour,stage,i) = imag(V .* (conj(Y) * conj(V)));
            end

            cost(i) = cost(i) + value(h)-dual(cons('VFA'));

            % write hourly storage values
            for hour = 1:param.H
                B_k(:, hour, stage,i) = max(dual(cons(sprintf('storage_%d', hour))),0);
            end

            % write hourly power generation
            for hour = 1:param.H
                P_g(:, hour, stage,i) = max(dual(cons(sprintf('generation_%d', hour))),0);
            end

            Q_g(:,:,stage,i) = -reshape(dual(cons('reactive')),m,param.H);

            % write curtailment values
            for hour = 1:param.H
                curtailment{i,1}.p_plus(:,hour,stage) = max(dual(cons(sprintf('curtail_p_plus_%d', hour))),0);
                curtailment{i,1}.p_minus(:,hour,stage) = max(dual(cons(sprintf('curtail_p_minus_%d', hour))),0);
                curtailment{i,1}.q_plus(:,hour,stage) = max(dual(cons(sprintf('curtail_q_plus_%d', hour))),0);
                curtailment{i,1}.q_minus(:,hour,stage) = max(dual(cons(sprintf('curtail_q_minus_%d', hour))),0);
            end
            % calculate resource state for the next stage
            new_sol(:,stage+1,i) = B_k(:,24,stage,i);
        end
    end
end

end
