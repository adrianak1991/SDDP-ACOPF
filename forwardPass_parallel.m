%%% This functions solves the OPF problems in forward pass for given
%%% scenarios where param.B_0*B_max is our initial value of storage.
%%% e0,e1,e2 contain parameters of cuts obtained until now.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[new_sol, obj_val, status, nodes, storage_hourly, curtailment_hourly] = forwardPass_parallel(data, param, num_scen, B_max, e0, e1, e2)

yalmip('clear')
[cons_bas, h_bas, dec] = defineProblem(data, param);

param.log_file = fopen(param.log_file_name, 'w');

% loading all network related data into the local workspace of the function
struct2vars(data.net);

% number of cuts presently in the VFA
num_cuts = length(e0{1,1});

% Resource states for all scenarios (initial levels at stage). 
% Since we also want to store the storage level at the end of planning, 
% i.e., after stage param.D, we need to store the initial level of stage param.D+1.
new_sol = zeros(k, param.D+1, num_scen); 
new_sol(:,1,:) = repmat(param.B_0*B_max, 1, 1, num_scen);

% hourly storage values (at the end of the respective hour)
storage_hourly = zeros(k, param.H, param.D, num_scen);

% hourly curtailment values
curtailment_hourly.p_plus = zeros(n, param.H, param.D, num_scen);
curtailment_hourly.p_minus = zeros(n, param.H, param.D, num_scen);
curtailment_hourly.q_plus = zeros(n, param.H, param.D, num_scen);
curtailment_hourly.q_minus = zeros(n, param.H, param.D, num_scen);

% collect values of objective function for all scenarios
obj_val = zeros(1,num_scen); 

% status information
status.errorCodes = zeros(param.D,param.H,num_scen); % Here we collect information if the condition for matrix A is satisfied in the given step
status.num_errors = zeros(param.D);
status.solvertime = 0; status.yalmiptime = 0;

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
nodes = ones(param.D+1,num_scen);
for stage = 1:param.D
    for i = 1:num_scen
        % simulate a node for the stage
        if stage == 1
            nodes(stage, i) = find(cumsum(data.prob{stage}(1,:)) >= rand, 1);
        else
            nodes(stage, i) = find(cumsum(data.prob{stage}(nodes(stage-1,i),:)) >= rand, 1);
        end
        
        % we add part of objective function involving the random parameter
        h = h_0 - dec.sigma(:,1)' * new_sol(:, stage, i);
        h = h + sum(dec.lambda .* squeeze(data.Pd_scen{stage}(nodes(i),:,:)),'all')+sum(dec.gamma .* squeeze(data.Qd_scen{stage}(nodes(i),:,:)),'all');
        
        % if there is already a VFA and its not the last stage (where there
        % is no VFA per definition)
        if num_cuts > 0 && stage < param.D
            h = h + e0{nodes(i),stage}' * nu;
            h = h + nu' * e1{nodes(i),stage} * B_max;
    
            % constraint with dual multiplier of last storage level
            cons = [cons_0, (e2{nodes(i), stage}' * nu + dec.sigma(:, param.H) + dec.kappa_u(:, param.H) >= 0):sprintf('storage_%d', param.H)];
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
                eig_val_vector = sort(diag(eig_val), 'ascend');
                
                W = dual(cons(sprintf('W_%d', hour)));
                [count, max_u, max_l, mean_dev, max_P, max_Q, ~, ~] = checkPrimal_power(data, param, eig_vec, double(dec.mu_u(:,hour)), double(dec.mu_l(:,hour)), W);
                
                if max_u > 1e-3 || max_l > 1e-3 || max_P > 1e-3 || max_Q > 1e-3
                  status.errorCodes(stage,hour,i) = 1;
                  status.num_errors(stage) = status.num_errors(stage) + 1;

                  fprintf(param.log_file, "F: %2d/%2d/%2d. acc: %5.2e \t (%f|%f|%f) (%f|%f|%f|%f|%f) \n", stage, i, hour, new_options.mosek.MSK_DPAR_INTPNT_CO_TOL_INFEAS, round(eig_val_vector(1),5), round(eig_val_vector(3),5), eig_val_vector(3)/eig_val_vector(2), count, max_u, max_l, mean_dev, max_P);
                end
            end

            % write hourly storage values
            for hour = 1:param.H
                storage_hourly(:, hour, stage, i) = max(dual(cons(sprintf('storage_%d', hour))),0);
            end
            
            % calculate resource state for the next stage (dual multiplier
            % of last constraint)
            new_sol(:,stage+1,i) = storage_hourly(:,param.H,stage,i);
            
            % write curtailment values
            for hour=1:param.H
                curtailment_hourly.p_plus(:,hour,stage,i) = max(dual(cons(sprintf('curtail_p_plus_%d', hour))),0);
                curtailment_hourly.p_minus(:,hour,stage,i) = max(dual(cons(sprintf('curtail_p_minus_%d', hour))),0);
                curtailment_hourly.q_plus(:,hour,stage,i) = max(dual(cons(sprintf('curtail_q_plus_%d', hour))),0);
                curtailment_hourly.q_minus(:,hour,stage,i) = max(dual(cons(sprintf('curtail_q_minus_%d', hour))),0);
            end
            
            % update objective value
            if num_cuts > 0 && stage < param.D
                obj_val(i) = obj_val(i) + value(h) - dual(cons('VFA'));
            else
                obj_val(i) = obj_val(i) + value(h);
            end
        end
    end
end

fclose(param.log_file);

end