%%% Here we implement modeling of load and simulation of given scenarios
%%% based on load residuals and change in scores from PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[scenarios, probability_transitions, N_final] = createLattice(data, param)

%% 1. READ DATA
%% Choose load data for specific area
if param.new_fit
    consumption = data.sets.histload(:,5);

    hydro = data.sets.histrenewables(:,12);
    PV = sum(data.sets.histrenewables(:,5:11),2) * param.VRES_mult;
    wind = data.sets.histrenewables(:,13) * param.VRES_mult;

    num_var = 4;
    data = [consumption PV hydro wind];

    %% 2. STOCHASTIC FACTORS
    %% Remove weekly and yearly seasonality

    % dates vector (data starts with a wednesday and is a leap year -> 2020)
    dates = datenum(2020,1,1,0,0,0):1/24:datenum(2020,12,31,23,0,0);
    data_residuals = zeros(size(data));

    % generate day dummies
    week_end = zeros(size(data,1), 1);
    day_type = { [2:5], [1,7] };
    I = ismember(weekday(dates), day_type{2});
    week_end(I) = 1;

    % generate hour dummies
    hours = zeros(size(data,1), 23);
    for h = 1:23
        hours(h:24:end,h) = 1;
    end

    % generate yearly seasonylity with two sine waves
    time = (1:size(data_residuals, 1))'; 
    yearly = [sin(2 * pi * time / 8784), cos(2 * pi * time / 8784)];

    % estiamte LASSO regression with 10 fold cross validation
    X = [ones(size(data,1), 1), week_end, hours, yearly];
    X2 = x2fx(X, 'quadratic');
    yearly_pattern = zeros(size(data, 1), num_var);
    residuals = zeros(size(data, 1), num_var);
    options = statset; options.UseParallel = true;
    beta = zeros(size(X2,2), num_var);
    a_0 = zeros(num_var, 1); a_1 = zeros(num_var, 1);
    for i=1:num_var
        [coeff, FitInfo] = lasso(X2, data(:,i), 'CV', 10, 'Options', options);
        beta(:,i) = coeff(:,FitInfo.IndexMinMSE);
        yearly_pattern(:,i) = max(X2 * beta(:,i) + FitInfo.Intercept(FitInfo.IndexMinMSE), 0);
        residuals(:,i) = data(:,i) - yearly_pattern(:,i);
        
        % quantile transform for residuals
        % Do a separate transformation for every hour of the day
        for h = 1:24
            sorted_residuals{i}{h} = unique(sort(residuals(h:24:end,i), 'ascend'));
            probabilities{i}{h} = ecdf(residuals(h:24:end,i));
            probabilities{i}{h} = probabilities{i}{h}(2:end);
            residuals(h:24:end,i) = interp1(sorted_residuals{i}{h}, probabilities{i}{h}, residuals(h:24:end,i));
        end
        
        fprintf('Variable: %d, R2: %f, MAPE: %f, MAD: %f\n', i, 1-var(yearly_pattern(:,i) - data(:,i)) / var(data(:,i)), mean(abs(yearly_pattern(:,i) - data(:,i)) ./ abs(yearly_pattern(:,i))), mean(abs(yearly_pattern(:,i) - data(:,i))));
    end
    
    % remove hydro modeling (not the best way to do it)
    residuals(:,3) = 0;
    yearly_pattern(:,3) = hydro;
    
    %% Analyze load residuals with PCA

    %% Conduct PCA
    residuals_daily = reshape(residuals', 24*num_var, numel(residuals) / (24*num_var))';

    [coeff, score, latent] = pca(residuals_daily, 'Centered', false);
    precision_levels = cumsum(latent') / sum(latent);
    number_components = find(precision_levels>=param.benchmark, 1);
    fprintf('Number of components in PCA: %d\n', number_components);

    % estimate VAR(1) model for score_changes using cross validated LASSO
    weights = latent(1:number_components);
    score = score(:,1:number_components) .* weights';
    lasso_residuals = zeros(size(score,1) - 1, size(score,2));
    lasso_coef = zeros(number_components, number_components);
    lasso_intercept = zeros(number_components, 1);
    for i = 1:number_components
        [B, FitInfo] = lasso(score(1:end-1, :), score(2:end, i), 'CV', 10);
        lasso_coef(i,:) = B(:, FitInfo.IndexMinMSE)';
        lasso_intercept(i) = FitInfo.Intercept(FitInfo.IndexMinMSE);
        lasso_residuals(:, i) = score(2:end, i)  -  score(1:end-1, :) * lasso_coef(i,:)' - lasso_intercept(i);
        fprintf('Lasso %d. R2: %f\n', i, 1-var(lasso_residuals(:,i))/var(score(:,i)));
    end
    
    save fit_data.mat;

    %% 3. SIMULATION
    S_0 = score(param.day_0, :)';
    sample_ids = randi(size(lasso_residuals, 1), param.num_scen * param.D, 1);
    sample_res = reshape(lasso_residuals(sample_ids,:)', [number_components, param.D, param.num_scen]);
    
    score_sample_lattice = zeros(number_components, param.D + 1, param.num_scen);
    score_sample_lattice(:, 1, :) = repmat(S_0, [1, 1, param.num_scen]);
    for i = 1:param.D
        score_sample_lattice(:, i+1, :) = lasso_coef * squeeze(score_sample_lattice(:, i, :)) + lasso_intercept + squeeze(sample_res(:,i,:));
    end
    
    %% 4. CREATING LATTICE
    addpath(genpath('quasar'));
    javaaddpath('quasar.jar');
    
    % change dimensions to make it fit to QUASAR (first number of scen, then
    % stages, then components)
    score_sample_lattice = permute(score_sample_lattice,[3,2,1]);
    
    % generate QUASAR time series sample
    sample_quasar = quasarTimeSeriesSample.createFromMultivariateSample(score_sample_lattice, cellstr(strcat('pc_',num2str((1:number_components)')))');
    
    % generate QUASAR lattice
    latBuilder = quasarLatticeBuilder.createFromSample(sample_quasar);
    latBuilder.numNodes(param.N);
    lattice_quasar = latBuilder.build;
    
    % save lattice to csv files
    lattice_quasar.storeCSV('myLattice');
else
    param_save = param;
    load fit_data.mat;
    param = param_save;
end

%% 5. READING DATA FROM CSV AND CREATING REQUIRED STRUCTURE
nodes_data = table2array(readtable('myLattice_nodes.csv','NumHeaderLines',1));
transitions_data = table2array(readtable('myLattice_transitions.csv','NumHeaderLines',1));

%% Create nodes of the lattice
nodes_numbers = cell(param.D, 1);
nodes_values = cell(param.D, 1);
N_final = param.N*ones(param.D,1);
scenarios = cell(param.D, 1);
for stage = 1:param.D
    nodes_in_stage = sum(nodes_data(:,1) == stage);
    nodes_numbers{stage} = nodes_data(nodes_data(:,1) == stage, 2);
    
    % first the components and then the nodes
    nodes_values{stage} = reshape(nodes_data(nodes_data(:,1) == stage, 4:end)', [number_components, nodes_in_stage]);
    
    % transform back to (24 x num_var) dimensional vector of quantile
    % transformed residuals
    data_residuals_matrix = (coeff(:,1:number_components)./weights')*nodes_values{stage};
    
    % round to [min(prob), 1] and then transform back to residuals using
    % interpolated inverse empirical CDF
    for i = 1:num_var
        for h = 1:24
            % transport data points back to the intervall [i/n,1], where i
            % is the number of occurances of the smallest residual for hour
            % h and variable i
            data_residuals_matrix((h-1) * num_var + i, :) = max(min(data_residuals_matrix((h-1) * num_var + i, :), 1), probabilities{i}{h}(1));
            
            % transform residuals back to natural range using inverse CDF
            data_residuals_matrix((h-1) * num_var + i, :) = interp1(probabilities{i}{h}, sorted_residuals{i}{h}, data_residuals_matrix((h-1) * num_var + i, :));
        end
    end
        
    % add seasonal trend
    new_data_residuals = reshape(data_residuals_matrix,[num_var, param.H, nodes_in_stage]);
    seasonal_correction = repmat(reshape(yearly_pattern((param.day_0+stage-1)*param.H+(1:param.H),:)',[num_var,param.H]),[1, 1, nodes_in_stage]);
    scenarios{stage} = new_data_residuals + seasonal_correction;
    scenarios{stage} = max(scenarios{stage},0);
    
    N_final(stage) = nodes_in_stage;
end


%% Create cell with probabilty transitions 
probability_transitions = cell(param.D,1);

first_transition = transitions_data(transitions_data(:,1)==0,2:3);
probability_transitions{1} = first_transition(first_transition(:,1)==nodes_numbers{1},2)';

for time = 1:(param.D-1)
    [p,q] = meshgrid(nodes_numbers{time+1}, nodes_numbers{time});
    required_values = [q(:) p(:)];
    [~,indices] = ismember(required_values,transitions_data(:,1:2),'rows');
    indices_reshape = reshape(indices, length(nodes_numbers{time}), length(nodes_numbers{time+1}));
    
    transition_values = zeros(length(nodes_numbers{time}), length(nodes_numbers{time+1}));
    transition_values(indices_reshape>0) = transitions_data(indices_reshape(indices_reshape>0),3);
    probability_transitions{time+1} = transition_values;
end

end