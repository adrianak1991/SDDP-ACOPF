% This function defines all decision variables, base constraints and base
% part of objective function for the primal problem on a node.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[cons_pbas,h_pbas,dec_p] = definePrimalProblem(data,param,B_max)

% transfer variables from struct to local workspace
struct2vars(data.net);
struct2vars(data.admittance);

% hourly generation
alpha_g = sdpvar(m,param.H,'full');
P_g = sdpvar(m,param.H,'full');
Q_g = sdpvar(m,param.H,'full');

% hourly storage values
B_k = sdpvar(k,param.H+1,'full');
C_plus = sdpvar(k,param.H,'full');
C_minus = sdpvar(k,param.H,'full');

% hourly curtailment values
DP_plus = sdpvar(n, param.H,'full');
DP_minus = sdpvar(n, param.H,'full');
DQ_plus = sdpvar(n, param.H,'full');
DQ_minus = sdpvar(n, param.H,'full');

% VF approximation
theta = sdpvar(1,1);

%% DEFINITION OF OBJECTIVE FUNCTION
h_pbas = sum(alpha_g,'all');
h_pbas = h_pbas + param.curtailment_cost_plus*sum(DP_plus+DQ_plus,'all') + param.curtailment_cost_minus*sum(DP_minus+DQ_minus,'all');
h_pbas = h_pbas + theta;

%% DEFINITION OF CONSTRAINTS
cons_pbas = [];

for time = 1:param.H
    % for cost function
    for j=1:(max_bp-1)
        cons_pbas = [cons_pbas, m_g(:,j).*(P_g(:,time)-a_g(:,j))+b_g(:,j)<=alpha_g(:,time)];
    end
    cons_pbas = [cons_pbas, P_min<=P_g(:,time), P_g(:,time)<=P_max];
    cons_pbas = [cons_pbas, Q_min<=Q_g(:,time), Q_g(:,time)<=Q_max];

    % storage
    cons_pbas = [cons_pbas, B_k(:,time+1)==B_k(:,time) + param.eta*C_plus(:,time)-C_minus(:,time)/param.eta];
    cons_pbas = [cons_pbas, B_k(:,time+1)<=B_max];
    cons_pbas = [cons_pbas, C_plus(:,time)<=B_max/(param.eta*4), C_minus(:,time)<=param.eta*B_max/4];
    cons_pbas = [cons_pbas, B_k(:,time+1)>=0, C_plus(:,time)>=0, C_minus(:,time)>=0];
    
    % curtailment
    cons_pbas = [cons_pbas,DP_plus(:,time)>=0, DP_minus(:,time)>=0, DQ_plus(:,time)>=0, DQ_minus(:,time)>=0];  
end

%% DEFINE A STRUCT WITH DECISION VARIABLES
dec_p.alpha_g = alpha_g; 
dec_p.P_g = P_g; 
dec_p.Q_g = Q_g;
dec_p.B_k = B_k;
dec_p.C_plus = C_plus;
dec_p.C_minus = C_minus;
dec_p.DP_plus = DP_plus;
dec_p.DP_minus = DP_minus; 
dec_p.DQ_plus = DQ_plus;
dec_p.DQ_minus = DQ_minus;
dec_p.theta = theta;

end