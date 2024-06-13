%%% Here we solve first stage optimization problem with investment decision
%%% on storage capacity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [storage_size, cost, objective, time] = firstStage(data, param, e0, e1)

% transfer variables from struct to local workspace
struct2vars(data.net);

storage_size = sdpvar(k, 1);
cons = [storage_size>=0, storage_size <= 1000];
h = sum(param.cost_storage.*storage_size);

if ~isempty(e0)
    theta0  = sdpvar;
    cons = [cons, e0 + e1 * storage_size <= theta0];
    h = h + theta0;
end

options = sdpsettings('solver','mosek','verbose', 0);
diagnostics = optimize(cons,h,options);

if diagnostics.problem~=0
    error('There was a problem in a forward pass of a first stage problem');
end

if ~isempty(e0)
    cost = value(h)-value(theta0);
else
    cost = value(h);
end

objective = value(h);
storage_size = value(storage_size);
time.yalmiptime = diagnostics.yalmiptime;
time.solvertime = diagnostics.solvertime;