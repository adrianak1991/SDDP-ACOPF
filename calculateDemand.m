%%% Here we implement the function which creates scenarios of demand
%%% for single buses from the data_scenarios for the whole area after
%%% subtracting renewable generation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Pd_scen,Qd_scen] = calculateDemand(data,param,data_scenarios)

n = data.net.n;
base = data.net.base;

bus_info = data.sets.bus_info;
N = param.N_final;

%% RENEWABLE GENERATION
bus_rg = data.net.bus_rg;
n_rg = length(bus_rg);

historical_PV = data.sets.histrenewables(:,5:11);
PV_ratios = sum(historical_PV,1)/sum(historical_PV,'all');

%% DEMAND
Pd_scen = cell(param.D,1);
Qd_scen = cell(param.D,1);

Pd_0 = bus_info(:,3);
Pd_0_area = sum(Pd_0);

Qd_0 = bus_info(:,4);

for stage=1:param.D
    renewables_scenarios = zeros(n_rg,param.H,N(stage));
    renewables_scenarios((n_rg-1):n_rg,:,:) = data_scenarios{stage}(3:4,:,:);% wind and hydro generation

    renewables_scenarios(1:(n_rg-2),:,:) = repmat(data_scenarios{stage}(2,:,:),n_rg-2,1,1).*repmat(PV_ratios',1,param.H,N(stage));

    Pd = repmat(data_scenarios{stage}(1,:,:),n,1,1).*repmat(Pd_0/Pd_0_area,1,param.H,N(stage));
    Pd(bus_rg,:,:) = Pd(bus_rg,:,:)-renewables_scenarios;
    Pd_scen{stage} = permute(Pd,[3,1,2])/base;

    Qd = repmat(data_scenarios{stage}(1,:,:),n,1,1).*repmat(Qd_0/Pd_0_area,1,param.H,N(stage));
    Qd_scen{stage} = permute(Qd,[3,1,2])/base;
end

end