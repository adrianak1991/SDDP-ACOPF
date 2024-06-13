%%% This function solves the primal problem for the provided active and
%%% reactive injections on a node and given scenarios. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cost,P_g,Q_g,B_k,curtailment]=optimalCostPrimal(data, param, num_scen, P, Q, B_max, nodes, e0, e1, e2)

% loading all network related data into the local workspace of the function
struct2vars(data.net);
struct2vars(data.admittance);

% define base components of optimization problem
[cons_pbas,h_pbas, dec_p] = definePrimalProblem(data, param, B_max);

% Resource states for all scenarios (initial levels at stage). since we also
% want to store the storage level of at the end of planning, i.e., after stage
% param.D, we need to store the initial level of stage param.D+1.
new_sol = zeros(k, param.D+1,num_scen); 
if size(B_max,2)==1
    new_sol(:,1,:) = repmat(param.B_0*B_max,1,num_scen);
else
    new_sol(:,1,:) = param.B_0*B_max;
end

cost =  sum(param.cost_storage.*B_max)*ones(num_scen,1);
P_g = zeros(m,param.H, param.D,num_scen);
Q_g = zeros(m,param.H, param.D,num_scen);
B_k = zeros(k,param.H,param.D,num_scen);

% hourly curtailment values
curtailment = cell(num_scen,1);
curtailment_cost_P = zeros(num_scen, 1);
curtailment_cost_Q = zeros(num_scen, 1);

%% PARAMETERS SPECIFIC FOR THIS PROBLEM
h = h_pbas;

%% RUNNING OPTIMIZATION PROBLEM FOR GIVEN SCENARIOS
for stage = 1:param.D
    for i=1:num_scen
        cons = cons_pbas;
        
        for hour=1:param.H
            % injections
            for bus = 1:n
                if(ismember(bus,storage_idx))
                    cons = [cons, P(bus,hour,stage,i) + data.Pd_scen{stage}(nodes(stage,i),bus,hour)-sum(dec_p.P_g(genidx(genidx(:,2)==bus,1),hour))+dec_p.C_plus(storage_idx==bus,hour)-dec_p.C_minus(storage_idx==bus,hour)-dec_p.DP_plus(bus,hour)+dec_p.DP_minus(bus,hour)==0];
                else
                    cons = [cons, P(bus,hour,stage,i) + data.Pd_scen{stage}(nodes(stage,i),bus,hour)-sum(dec_p.P_g(genidx(genidx(:,2)==bus,1),hour))-dec_p.DP_plus(bus,hour)+dec_p.DP_minus(bus,hour)==0];
                end
                
                cons = [cons, Q(bus,hour,stage,i) + data.Qd_scen{stage}(nodes(stage,i),bus,hour)-sum(dec_p.Q_g(genidx(genidx(:,2)==bus,1),hour))-dec_p.DQ_plus(bus,hour)+dec_p.DQ_minus(bus,hour)==0];  
            end
        end
        
        new_sol(new_sol(:,stage,i) < 1e-5,stage,i) = 0;
        cons = [cons, dec_p.B_k(:,1)==new_sol(:,stage,i)]; 
        if stage<param.D
            cons = [cons, e0{nodes(stage,i),stage} + e1{nodes(stage,i),stage} * B_max + e2{nodes(stage,i),stage} * dec_p.B_k(:,param.H+1)<=dec_p.theta];
        else
            cons = [cons, dec_p.B_k(:,param.H+1)==param.B_0*B_max, 0<=dec_p.theta];
        end
        
        options  =  sdpsettings('solver','mosek','verbose',0);
        diagnostics = optimize(cons, h, options);
        
        if(diagnostics.problem ~= 0)
            error(stage + " Primal: Solver error: " + diagnostics.info);
        else
            cost(i) = cost(i) + value(h)-value(dec_p.theta);

            % write hourly storage values
            for hour = 1:param.H
                B_k(:, hour, stage,i) = max(value(dec_p.B_k(:,hour+1)),0);
            end

            % write hourly power generation
            P_g(:,:, stage,i) = max(value(dec_p.P_g),0);
            Q_g(:,:,stage,i) =  value(dec_p.Q_g);

            % write curtailment values
            curtailment{i,1}.p_plus(:,:,stage) = max(value(dec_p.DP_plus),0);
            curtailment{i,1}.p_minus(:,:,stage) = max(value(dec_p.DP_minus),0);
            curtailment{i,1}.q_plus(:,:,stage) = max(value(dec_p.DQ_plus),0);
            curtailment{i,1}.q_minus(:,:,stage) = max(value(dec_p.DQ_minus),0);
            
            % Calculate resource state for the next stage (dual multiplier
            % of last constraint)
            new_sol(:,stage+1,i) = B_k(:, param.H, stage,i);

            % curtailment cost
            curtailment_cost_P(i) = curtailment_cost_P(i) + double(param.curtailment_cost_plus*sum(dec_p.DP_plus,'all')+param.curtailment_cost_minus*sum(dec_p.DP_minus,'all'));
            curtailment_cost_Q(i) = curtailment_cost_Q(i) + double(param.curtailment_cost_plus*sum(dec_p.DQ_plus,'all')+param.curtailment_cost_minus*sum(dec_p.DQ_minus,'all'));
        end
    end
end

end
