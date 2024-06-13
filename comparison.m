%% Finding a solution for certain scenarios
param.seed = 5342; rng(param.seed);
num_scen = 100; % how many scenarios we generate

assignment_scen = cell(param.num_workers, 1);
s_scen = floor(num_scen/param.num_workers);
for w = 1:param.num_workers
    assignment_scen{w} = s_scen*(w-1)+1:s_scen*w;

    if w <= mod(num_scen, param.num_workers)
        assignment_scen{w} = [assignment_scen{w}, param.num_workers*s_scen+w];
    end
end

% We generate some scenarios and compare costs obtained by solving dual
% problem with optimal cost from primal problem for given voltages
param.EV_4D = true;
spmd
    [cost, W, P, Q, P_g, Q_g, B_k, curtailment, B_max, nodes, status] = findSolution(data, param, length(assignment_scen{labindex,1}), e0_0, e1_0, e0, e1, e2);
    [cost_p,P_g_p,Q_g_p,B_k_p,curtailment_p] = optimalCostPrimal(data, param, length(assignment_scen{labindex,1}), P, Q, B_max,nodes, e0, e1, e2);
end
cost = vertcat(cost{:});
cost_p = vertcat(cost_p{:});
W = cat(3,W{:});
P_g = cat(4,P_g{:});
Q_g = cat(4,Q_g{:});
B_k = cat(4,B_k{:});
B_max = [B_max{:}];
nodes = [nodes{:}];
curtailment = vertcat(curtailment{:});
curtailment_p = vertcat(curtailment_p{:});

% Plot storage decisions over time and the associated residual demand levels at the node
storage = 8;
storage_level = reshape(squeeze(B_k(storage, :, :, :)), param.H * param.D, num_scen);
figure; subplot(2,1,1); fanChart(1:param.D*param.H, storage_level*100); title('Storage Levels');xlim([0,168]);
xlabel("Time (h)"); ylabel("Storage level (MWh)");
set(gca,'FontSize',16)

node = data.net.storage_idx(storage);
residual_demand = zeros(param.H * param.D, num_scen);
for scen = 1:num_scen
    for stage = 1:param.D
        residual_demand((stage-1)*24+1:stage*24, scen) = squeeze(data.Pd_scen{stage}(nodes(stage, scen), node, :));
    end
end
subplot(2,1,2); fanChart(1:param.D*param.H, residual_demand*100); title('Residual Active Power Demand');xlim([0,168]);
xlabel("Time (h)"); ylabel("Demand (MW)");
set(gca,'FontSize',16)


% Plot curtailment values (sum up over all buses) for dual and primal solution
curtailment_sum = zeros(param.H*param.D,num_scen,2); %plus, minus
curtailment_p_sum = zeros(param.H*param.D,num_scen,2);

for i=1:num_scen
   curtailment_sum(:,i,1) = reshape(sum(curtailment{i,1}.p_plus,1) + sum(curtailment{i,1}.q_plus,1),param.H*param.D,1,1);
   curtailment_sum(:,i,2) = reshape(sum(curtailment{i,1}.p_minus,1) + sum(curtailment{i,1}.q_minus,1),param.H*param.D,1,1);
   curtailment_p_sum(:,i,1) = reshape(sum(curtailment_p{i,1}.p_plus,1) + sum(curtailment_p{i,1}.q_plus,1),param.H*param.D,1,1);
   curtailment_p_sum(:,i,2) = reshape(sum(curtailment_p{i,1}.p_minus,1)+ sum(curtailment_p{i,1}.q_minus,1),param.H*param.D,1,1);
end

for i = 1:num_scen
    fprintf('DUAL curtailment. Scenario %d, P_plus: %.2f, P_minus: %.2f, Q_plus: : %.2f, Q_minus: : %.2f\n', i, sum(curtailment{i}.p_plus,'all') * param.curtailment_cost_plus, sum(curtailment{i}.p_minus,'all') * param.curtailment_cost_minus, sum(curtailment{i}.q_plus,'all') * param.curtailment_cost_plus, sum(curtailment{i}.q_minus,'all') * param.curtailment_cost_minus);
    fprintf('PRIMAL curtailment. Scenario %d, P_plus: %.2f, P_minus: %.2f, Q_plus: : %.2f, Q_minus: : %.2f\n', i, sum(curtailment_p{i}.p_plus,'all') * param.curtailment_cost_plus, sum(curtailment_p{i}.p_minus,'all') * param.curtailment_cost_minus, sum(curtailment_p{i}.q_plus,'all') * param.curtailment_cost_plus, sum(curtailment_p{i}.q_minus,'all') * param.curtailment_cost_minus);
end
disp("Average deviations between optimal value of dual problem and primal for given voltages: " +sum(cost_p-cost)/sum(cost_p)); %average deviations
toc
