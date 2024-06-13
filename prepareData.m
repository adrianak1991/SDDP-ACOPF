%%% Here we implement the function which prepare all data about the network
%%% required to solve the problem and put together in one data structure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[data_struct] = prepareData(area)

%% GETTING AND TRANSFORMATION OF DATA
%% - getting raw data
case_data = case_RTS_GMLC;
renewables_data = importdata('DAY_AHEAD_renewables_area1.csv',';',1);
load_data = importdata('DAY_AHEAD_regional_Load.csv',',',1); 
gen_data = readtable('gen.csv'); % information about generators (extended version from github)
gen_data = gen_data(:,[2 7]); % we extract number of bus and fuel type of generator

%% - creating set of all types of generators in form of vector and adding 
% information about type to the main set (represented now by number)
gen_types = unique(gen_data(:,2));

[wasfound,idx] = ismember(gen_data(:,2),gen_types);
case_data.gen = [case_data.gen idx];
gen_types = string(gen_types.Variables);

%% - reducing network to the given area (update of all sets from case_data)
area_lb = case_data.areas(area,2);
area_ub = case_data.areas(area+1,2);

case_data.bus = case_data.bus(area_lb<=case_data.bus(:,1) & case_data.bus(:,1)<area_ub,:);

case_data.branch = case_data.branch(area_lb<=case_data.branch(:,1) & case_data.branch(:,1)<area_ub,:);
case_data.branch = case_data.branch(area_lb<=case_data.branch(:,2) & case_data.branch(:,2)<area_ub,:);

i = area_lb<=case_data.gen(:,1) & case_data.gen(:,1)<area_ub;
case_data.gen = case_data.gen(i,:);
case_data.gencost = case_data.gencost(i,:);

load_data = load_data.data(:,[1:4 4+area]);

%% MODIFICATION OF SET OF GENERATORS
%% - removal of generators out of service
i = find(case_data.gen(:,8)==0);
case_data.gen(i,:) = [];
case_data.gencost(i,:) = [];

%% - removal of generators with zero cost (renewable power plants are included in residual demand)
i = find(case_data.gencost(:,6)==0 & case_data.gencost(:,8)==0);
case_data.gen(i,:) = [];
case_data.gencost(i,:) = [];

%% - removal of coal power plants
% remove 657 MW of coal capacity (59% of total capacity)
I = find(case_data.gen(:,end)==find(gen_types=='Oil'));
case_data.gen(I,:) = [];
case_data.gencost(I,:) = [];

J = find(case_data.gen(:, 1) == 101 & case_data.gen(:,end)==find(gen_types=='Coal'), 1);
J = [J, find(case_data.gen(:, 1) == 102 & case_data.gen(:,end)==find(gen_types=='Coal'), 1)];
J = [J, find(case_data.gen(:, 1) == 123 & case_data.gen(:, 2) == 350 & case_data.gen(:,end)==find(gen_types=='Coal'), 1)];
J = [J, find(case_data.gen(:, 1) == 115 & case_data.gen(:,end)==find(gen_types=='Coal'), 1)];
J = [J, find(case_data.gen(:, 1) == 116 & case_data.gen(:,end)==find(gen_types=='Coal'), 1)];
case_data.gen(J,:) = [];
case_data.gencost(J,:) = [];

%% - setting minimum generation to zero
case_data.gen(:,10) = zeros(length(case_data.gen(:,1)),1);

%% - introduction of bus numbers and generators numbers in all data sets
n = length(case_data.bus(:,1));     % number of all buses
m = length(case_data.gen(:,1));     % number of generators
l = length(case_data.branch(:,1));  % number of lines including directions

%% - introduction of new numbers for buses
num_trans = [(1:n)' case_data.bus(:,1)]; % new number, old number
case_data.bus(:,1) = (1:n)';

for i=1:m
    case_data.gen(i,1) = find(num_trans(:,2)==case_data.gen(i,1));
end
case_data.gen = [(1:m)' case_data.gen];
case_data.gencost = [(1:m)' case_data.gencost];

for i=1:l
    case_data.branch(i,1) = find(num_trans(:,2)==case_data.branch(i,1));
    case_data.branch(i,2) = find(num_trans(:,2)==case_data.branch(i,2));
end

%% - modification of bus numbers in data set with renewables generation
headers = str2double(renewables_data.colheaders(5:end));

for i=1:length(headers)
    renewables_data.colheaders(4+i) = {find(num_trans(:,2)==headers(i))};
end

%% PREPARATION OF DATA STRUCTURE
%% - saving data sets we use
data_struct.sets.histrenewables = renewables_data.data;
data_struct.sets.histload = load_data;
data_struct.sets.bus_info = case_data.bus;
data_struct.sets.branch_info = case_data.branch;
data_struct.sets.gen_info = case_data.gen;
data_struct.sets.gen_cost = case_data.gencost;

%% - saving basic parametrs of network
data_struct.net.base = case_data.baseMVA;
data_struct.net.n = length(case_data.bus(:,1));
data_struct.net.n_0 = case_data.bus(case_data.bus(:,2)==3,1);
data_struct.net.m = length(case_data.gen(:,1));
data_struct.net.l = length(case_data.branch(:,1));
data_struct.net.L = case_data.branch(:,1:2);
data_struct.net.genidx = case_data.gen(:,1:2);

data_struct.net.bus_rg = cell2mat(renewables_data.colheaders(5:end)); % buses with renewable generation
data_struct.net.n_rg = length(data_struct.net.bus_rg);
data_struct.net.gen_types = gen_types;

%% - parameters of generators and buses
data_struct.net.P_min = case_data.gen(:,11)/case_data.baseMVA;
data_struct.net.P_max = case_data.gen(:,10)/case_data.baseMVA; % upper bound on real power generation

data_struct.net.Q_min = case_data.gen(:,6)/case_data.baseMVA; % lower bound on reactive power generation
data_struct.net.Q_max = case_data.gen(:,5)/case_data.baseMVA; % upper bound on reactive power generation

data_struct.net.V_min = case_data.bus(:,13); % lower bound on voltage magnitude
data_struct.net.V_max = case_data.bus(:,12); % upper bound on voltage magnitude

%% - parameters of branches
% upper bound on apparent power transferred through the line
data_struct.net.Sline_max = repmat(case_data.branch(:,6),2,1)/case_data.baseMVA;
data_struct.net.Pline_max = repmat(case_data.branch(:,6),2,1)/case_data.baseMVA;

%% - parameters of cost function
% here we define cost function for generators which have piecewise linear cost function
cost_function_pw = case_data.gencost(case_data.gencost(:,2)==1,:);
data_struct.net.genidx_pw = cost_function_pw(:,1);  

% preparation data for piecewise linear cost function
ng_pw = length(cost_function_pw(:,1)); % number of generators with piecewise linear cost function
n_bp = cost_function_pw(:,5); % number of break points
max_bp = max(n_bp);

a_g = cost_function_pw(:,5+(1:2:2*max_bp))/case_data.baseMVA;
b_g = cost_function_pw(:,5+(2:2:2*max_bp));
m_g = (b_g(:,2:max_bp)-b_g(:,1:(max_bp-1)))./(a_g(:,2:max_bp)-a_g(:,1:(max_bp-1)));

a_g2 = [zeros(ng_pw,1) a_g(:,2:(max_bp-1))];
b_g2 = zeros(ng_pw,max_bp-1);
for i=2:(max_bp-1)
    b_g2(:,i) = m_g(:,i-1).*(a_g2(:,i)-a_g2(:,i-1))+b_g2(:,i-1);
end

data_struct.net.cost_function_pw = cost_function_pw;
data_struct.net.ng_pw = ng_pw;
data_struct.net.n_bp = n_bp;
data_struct.net.max_bp = max_bp;
data_struct.net.m_g = m_g;
data_struct.net.a_g = a_g2;
data_struct.net.b_g = b_g2;

%% - parameters of storage
storage_idx = unique(data_struct.net.bus_rg)';  % indexes of buses with storage
k = length(storage_idx);     % number of storages

data_struct.net.k = k;
data_struct.net.storage_idx = storage_idx;

%% - creating admittance matrix
data_struct.admittance = createAdmittanceMatrix(data_struct);

end
