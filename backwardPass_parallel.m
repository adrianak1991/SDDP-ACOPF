% This function solves one step of backward pass for given set of scenarios
% by adding the changing part of problem to the base components.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[e0_val,e1_val,e2_val,status] = backwardPass_parallel(data, param, scenarios, stage, B_max, sol, e0, e1, e2)

yalmip('clear')
[cons_bas, h_bas, dec] = defineProblem(data, param);

param.log_file = fopen(param.log_file_name, 'w');

% loading a file with all parameters required to solve the problem (i.e. all limits).
struct2vars(data.net);

% number of cuts to be computed
num_scen = size(scenarios,1);

% number of cuts presently in the VFA in the given stage
num_cuts = length(e0{1,stage});

%% Here we initialize values to be returned.
e0_val = zeros(1,num_scen);
e1_val = zeros(k,num_scen);
e2_val = zeros(k,num_scen);
status.errorCodes = zeros(param.H,num_scen);
status.num_errors = 0;
status.yalmiptime = 0; status.solvertime = 0;

%% PARAMETERS SPECIFIC FOR THIS PROBLEM
h0 = h_bas - sum((dec.kappa_u + dec.beta_u_minus * param.eta/4 + dec.beta_u_plus / (4*param.eta))' * B_max);
h0 = h0 - dec.sigma(:,1)' * sol';

% only important if there is a VFA already
if num_cuts > 0 && stage < param.D
    % Lagrange multiplier associated with value function approximation
    % constraint.
    nu = sdpvar(num_cuts, 1);
    cons0 = [cons_bas, nu>=0, (sum(nu)==1):'VFA'];
else
    cons0 = cons_bas;
end

if stage==param.D
    tau = sdpvar(k,1,'full');
    
    h0 = h0 - tau'*(param.B_0*B_max);
end

%% CHANGING PART OF THE PROBLEM
% Solving the problem for every scenario at the given stage
for i=1:num_scen
    h = h0 + sum(dec.lambda .* squeeze(data.Pd_scen{stage}(scenarios(i),:,:)), 'all') + sum(dec.gamma .* squeeze(data.Qd_scen{stage}(scenarios(i),:,:)),'all');
    
    % distinguish between the first and subsequent iterations
    if num_cuts > 0 && stage < param.D
        h = h + e0{scenarios(i),stage}' * nu;
        h = h + nu' * e1{scenarios(i),stage} * B_max;
        cons = [cons0, e2{scenarios(i), stage}' * nu + dec.sigma(:, param.H) + dec.kappa_u(:, param.H) >= 0];
    else
        if stage==param.D
            cons = [cons0, dec.sigma(:, param.H) + dec.kappa_u(:, param.H) + tau >= 0];
        else
            cons = [cons0, dec.sigma(:, param.H) + dec.kappa_u(:, param.H) >= 0];
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

    if(diagnostics.problem~=0)
        error("B: Solver error: " + diagnostics.info);
    else
        % We have to check if required condition is satisfied (number of
        % zero-eigenvalues is equal to 2).
        for hour = 1:param.H
            A_sol = value(dec.matA(:,:,hour));
            [eig_vec, eig_val] = eig(A_sol);
            eig_val_vector = sort(diag(eig_val), 'ascend');
            
            W = dual(cons(sprintf('W_%d', hour)));
            [count, max_u, max_l, mean_dev, max_P, max_Q, ~, ~] = checkPrimal_power(data, param, eig_vec, double(dec.mu_u(:,hour)), double(dec.mu_l(:,hour)), W);

            if max_u > 1e-3 || max_l > 1e-3 || max_P > 1e-3 || max_Q > 1e-3
              status.errorCodes(hour,i) = 1;
              status.num_errors = status.num_errors + 1;

              fprintf(param.log_file, "B: %2d/%2d/%2d. acc: %5.2e \t (%f|%f|%f) (%f|%f|%f|%f|%f) \n", stage, i, hour, new_options.mosek.MSK_DPAR_INTPNT_CO_TOL_INFEAS, round(eig_val_vector(1),5), round(eig_val_vector(3),5), eig_val_vector(3)/eig_val_vector(2), count, max_u, max_l, mean_dev, max_P);
            end            
        end
           
        % calculate slope and intercept of cut for the scenario
        e2_val(:,i) = -value(dec.sigma(:,1));
        if num_cuts > 0 && stage < param.D
            e1_val(:,i) =  e1{scenarios(i),stage}' * value(nu) - sum(value(dec.kappa_u) + value(dec.beta_u_minus)*param.eta/4 + value(dec.beta_u_plus) / (4*param.eta), 2);
        else
            if stage==param.D
                e1_val(:,i) =  -sum(value(dec.kappa_u) + value(dec.beta_u_minus)*param.eta/4 + value(dec.beta_u_plus) / (4*param.eta), 2)-value(tau)*param.B_0;
            else
                e1_val(:,i) =  -sum(value(dec.kappa_u) + value(dec.beta_u_minus)*param.eta/4 + value(dec.beta_u_plus) / (4*param.eta), 2);
            end
        end
        e0_val(i) =  value(h) - e1_val(:,i)' * B_max - e2_val(:,i)'*sol';
    end
end

fclose(param.log_file);

end