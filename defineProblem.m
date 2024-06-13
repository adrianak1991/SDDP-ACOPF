% This function defines all decision variables, base constraints and base
% part of objective function for the dual problem on a node. The idea is to 
% reuse these variables and constraints to avoid repeated definition during
% the running of the SDDP algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[cons_bas,h_bas, dec] = defineProblem(data, param)

% transfer variables from struct to local workspace
struct2vars(data.net);
struct2vars(data.admittance);

% basic constraints 
cons_bas = [];

%% DECISION VARIABLES
% Lagrange multiplier associated with piecewise linear cost function.
zeta = sdpvar(ng_pw,max_bp-1,param.H,'full');
cons_bas = [cons_bas, zeta>=0];

% Lagrange multipliers associated with upper bound of active power for a
% generator.
lambda_u = sdpvar(m, param.H, 'full');
cons_bas = [cons_bas, lambda_u >= 0];

% Lagrange multiplier associated with definition of active power on a node
lambda = sdpvar(n, param.H, 'full');

% Lagrange multipliers associated with lower and upper bound of reactive 
% power provided by generators.
gamma_l = sdpvar(m, param.H, 'full');
gamma_u = sdpvar(m, param.H, 'full');
cons_bas = [cons_bas, gamma_u>=0, gamma_l>=0];

% Lagrange multiplier associated with definition of reactive power at the
% buses of the network.
gamma = sdpvar(n, param.H, 'full');

% Lagrange multipliers associated with lower and upper bounds for voltage.
mu_l = sdpvar(n, param.H, 'full');
mu_u = sdpvar(n, param.H, 'full');
cons_bas = [cons_bas, mu_u>=0, mu_l>=0];

% Lagrange multiplier associated with matrix representing upper bound of
% apparent power flow.
r_lm = sdpvar(2*l, 5, param.H, 'full');

% Lagrange multipliers associated with defition of storage level.
sigma = sdpvar(k, param.H, 'full');

% Lagrange multipliers associated with upper bound of storage level.
kappa_u = sdpvar(k, param.H, 'full');
cons_bas = [cons_bas, kappa_u>=0];

% Lagrange multipliers associated with upper bound of real power 
% charged and discharged into storage.
beta_u_plus = sdpvar(k, param.H, 'full');
beta_u_minus = sdpvar(k, param.H, 'full');
cons_bas = [cons_bas, beta_u_plus>=0, beta_u_minus>=0];

%% DEFINE BASE OBJECTIVE FUNCTION AND CONSTRAINTS
matA = sdpvar(2*n, 2*n, param.H, 'full'); 
 
%% OBJECTIVE FUNCTION
% Define part of the objective function that doesn't depend on the resource
% state, random variables realizations or the value function approximations
h_bas = sum(zeta.*repmat(b_g-m_g.*a_g, 1, 1, param.H), 'all');
h_bas = h_bas + sum(-lambda_u' * P_max);
h_bas = h_bas + sum(gamma_l' * Q_min - gamma_u' * Q_max);
h_bas = h_bas + sum(mu_l' * (V_min.^2) - mu_u' * (V_max.^2));
h_bas = h_bas - (Sline_max.^2)'*sum(r_lm(:,1,:),3)-sum(r_lm(:,4:5,:),'all');
    
for time = 1:param.H
    % Define matrix A
    A = sum(repmat(lambda(:,time),1,2*n,2*n).*bY_k,1);
    A = A + sum(repmat(gamma(:,time),1,2*n,2*n).*hbY_k,1);
    A = A + sum(repmat((mu_u(:,time)-mu_l(:,time)),1,2*n,2*n).*M_k,1);

    A = A + sum(repmat(2*r_lm(:,2,time),1,2*n,2*n).*bY_lm,1);
    A = A + sum(repmat(2*r_lm(:,3,time),1,2*n,2*n).*hbY_lm,1);

    A = squeeze(A);

    % Define constraints
    % store matrices A for all hours and impose psd constraint
    cons_bas = [cons_bas, (A>=0):sprintf('W_%d', time), matA(:,:,time)==A]; 
    cons_bas = [cons_bas, (sum(zeta(:,:,time).*m_g,2)+lambda_u(genidx_pw,time)-lambda(genidx(genidx_pw,2),time)>=0):sprintf('generation_%d', time)];

    % psd constraints for dual multipliers of active power on line constraints
    for j=1:2*l
        cons_bas = [cons_bas, [r_lm(j,1:3,time); r_lm(j,2,time) r_lm(j, 4,time) 0; r_lm(j,3,time) 0 r_lm(j,5,time)]>=0];
    end
    
    % storage
    if time < param.H
        cons_bas = [cons_bas, (sigma(:,time)-sigma(:,time+1)+kappa_u(:,time)>=0):sprintf('storage_%d', time)];
    end
    
    % curtailment
    cons_bas = [cons_bas, (param.curtailment_cost_plus - lambda(:,time)>=0):sprintf('curtail_p_plus_%d', time)];
    cons_bas = [cons_bas, (param.curtailment_cost_minus + lambda(:,time)>=0):sprintf('curtail_p_minus_%d', time)];
    cons_bas = [cons_bas, (param.curtailment_cost_plus - gamma(:,time)>=0):sprintf('curtail_q_plus_%d', time)];
    cons_bas = [cons_bas, (param.curtailment_cost_minus + gamma(:,time)>=0):sprintf('curtail_q_minus_%d', time)];
 
end

% zetas have to sum up to one for all generators (m) and all hours
cons_bas = [cons_bas, sum(zeta,2)==ones(m,1,param.H)]; 
% reactive power
cons_bas = [cons_bas, (gamma_u-gamma_l-gamma(genidx(:,2),:)==0):'reactive'];
% storage injection
cons_bas = [cons_bas, lambda(storage_idx,:) - sigma * param.eta + beta_u_plus>=0];
% storage withdrawal
cons_bas = [cons_bas, -lambda(storage_idx,:) + sigma / param.eta + beta_u_minus>=0];


%% DEFINE A STRUCT WITH DECISION VARIABLES
dec.zeta = zeta;
dec.lambda_u = lambda_u;
dec.lambda = lambda;
dec.gamma_l = gamma_l;
dec.gamma_u = gamma_u;
dec.gamma = gamma;
dec.mu_l = mu_l;
dec.mu_u = mu_u;
dec.r_lm = r_lm;
dec.sigma = sigma;
dec.kappa_u = kappa_u;
dec.beta_u_plus = beta_u_plus;
dec.beta_u_minus = beta_u_minus;
dec.matA = matA;