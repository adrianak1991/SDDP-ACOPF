% Calculate storage and plant operations given a power flow from a matrix A
function [count, max_u, max_l, mean_dev, max_P, max_Q, V, curtailment_cost] = checkPrimal_power(data, param, eig_vec, mu_u, mu_l, W)
    struct2vars(data.net);
    struct2vars(data.admittance);
    curtailment_cost = 0;
    
    % Compute power balance on nodes from matrix W
    P_D = zeros(n,1); Q_D = zeros(n,1);
    for i = 1:n
        P_D(i) = trace(squeeze(bY_k(i,:,:))' * W);
        Q_D(i) = trace(squeeze(hbY_k(i,:,:))' * W);
    end

    %%%% 1. Reconstruct voltages and power
    F1 = eig_vec(1:n,1);
    F2 = eig_vec(n+1:end, 1);

    % find all binding voltage constraints
    J_u = find(mu_u > 10^(-5)); J_l = find(mu_l > 10^(-5)); max_error = inf;
    if length(J_u) + length(J_l) == 0
        zeta = [F1(n_0) -F2(n_0); F2(n_0) F1(n_0)] \ [1; 0]; % that is not a good alternative

        V = complex(zeta(1) * F1 - zeta(2) * F2, zeta(1) * F2 + zeta(2) * F1);
        I = abs(V) < V_min | abs(V) > V_max;
        max_l = max(max(V_min - abs(V)),0); max_u = max(max(abs(V) - V_max), 0);
        count = sum(I); mean_dev = mean(max(V_min - abs(V), 0) + max(abs(V) - V_max, 0));
        max_error = max(max_l, max_u);

        % check errors in the power balance of every node
        P_V = real(V .* (conj(Y) * conj(V))); % active power as calculated from the voltages
        Q_V = imag(V .* (conj(Y) * conj(V))); % reactive power as calculated from the voltages
        max_P = max(abs(P_V - P_D));
        max_Q = max(abs(Q_V - Q_D));
        
        curtailment_cost = sum((max(P_V - P_D,0) + max(Q_V - Q_D,0))) * param.curtailment_cost_plus - sum((min(P_V - P_D,0) + min(Q_V - Q_D,0))) * param.curtailment_cost_minus;
        fprintf(param.log_file, 'WARNING: No binding voltage constraint. max_l: %g, max_u: %g, max_error: %g\n', max(mu_l), max(mu_u), max_error);
    else
        for i = 1:length(J_u)
            zeta(2) = sqrt( V_max(J_u(i))^2 / ((F1(n_0)^2/F2(n_0)^2 + 1)*(F1(J_u(i))^2+F2(J_u(i))^2))); % minus of that would also be possible
            zeta(1) = -F1(n_0)/F2(n_0) * zeta(2);

            V_test = complex(zeta(1) * F1 - zeta(2) * F2, zeta(1) * F2 + zeta(2) * F1);
            I = abs(V_test) < V_min | abs(V_test) > V_max;

            % compute errors in the upper and lower bounds of voltages in
            % all nodes
            max_l_test = max(max(V_min - abs(V_test)),0); max_u_test = max(max(abs(V_test) - V_max), 0);

            % check errors in the power balance of every node
            P_V = real(V_test .* (conj(Y) * conj(V_test))); % active power as calculated from the voltages
            Q_V = imag(V_test .* (conj(Y) * conj(V_test))); % reactive power as calculated from the voltages
            max_P_test = max(abs(P_V - P_D));
            max_Q_test = max(abs(Q_V - Q_D));

            if max([max_l_test, max_u_test, max_P_test, max_Q_test]) < max_error
                max_error = max([max_l_test, max_u_test, max_P_test]);
                max_u = max_u_test;
                max_l = max_l_test;
                max_P = max_P_test;
                max_Q = max_Q_test;
                V = V_test;

                count = sum(I); mean_dev = mean(max(V_min - abs(V), 0) + max(abs(V) - V_max, 0));
            end
        end

        for i = 1:length(J_l)
            zeta(2) = sqrt( V_min(J_l(i))^2 / ((F1(n_0)^2/F2(n_0)^2 + 1)*(F1(J_l(i))^2+F2(J_l(i))^2))); % minus of that would also be possible
            zeta(1) = -F1(n_0)/F2(n_0) * zeta(2);

            V_test = complex(zeta(1) * F1 - zeta(2) * F2, zeta(1) * F2 + zeta(2) * F1);
            I = abs(V_test) < V_min | abs(V_test) > V_max;

            % compute errors in the upper and lower bounds of voltages in
            % all nodes
            max_l_test = max(max(V_min - abs(V_test)),0); max_u_test = max(max(abs(V_test) - V_max), 0);

            % check errors in the power balance of every node
            P_V = real(V_test .* (conj(Y) * conj(V_test))); % active power as calculated from the voltages
            Q_V = imag(V_test .* (conj(Y) * conj(V_test))); % reactive power as calculated from the voltages
            max_P_test = max(abs(P_V - P_D));
            max_Q_test = max(abs(Q_V - Q_D));

            if max([max_l_test, max_u_test, max_P_test, max_Q_test]) < max_error
                max_error = max([max_l_test, max_u_test, max_P_test, max_Q_test]);
                max_u = max_u_test;
                max_l = max_l_test;
                max_P = max_P_test;
                max_Q = max_Q_test;
                V = V_test;

                count = sum(I); mean_dev = mean(max(V_min - abs(V), 0) + max(abs(V) - V_max, 0));
            end
        end
        
        curtailment_cost = sum((max(real(V .* (conj(Y) * conj(V))) - P_D,0) + max(imag(V .* (conj(Y) * conj(V))) - Q_D,0))) * param.curtailment_cost_plus - sum((min(real(V .* (conj(Y) * conj(V))) - P_D,0) + min(imag(V .* (conj(Y) * conj(V))) - Q_D,0))) * param.curtailment_cost_minus;
    end

    % if there are more than 2 non-zero eigenvalues and the maxmimal deviations
    % are too large, try to fix the problem by using 4 conditions
    if param.EV_4D && max_error > 1e-3
        error_save = max_error;
        V_save = V;
        V_0 = [real(V); imag(V)];

        % find voltages that fit power at nodes best (respecting upper and
        % lower limits)
        opt = optimoptions('fmincon');
        opt.MaxIterations = 10000; opt.MaxFunctionEvaluations = 1e6; opt.Display = 'off';
        opt.Algorithm = 'sqp';
        [V, error] = fmincon(@(V) mean_deviation(V, P_D, Q_D, Y), V_0, [], [], [], [], [], [], @(V) voltage_constraints(V, V_max, V_min), opt);

        V_c = complex(V(1:length(V)/2), V(length(V)/2+1:end));
        P_V = real(V_c .* (conj(Y) * conj(V_c))); % active power as calculated from the voltages
        Q_V = imag(V_c .* (conj(Y) * conj(V_c))); % reactive power as calculated from the voltages
        max_P_test = max(abs(P_V - P_D));
        max_Q_test = max(abs(Q_V - Q_D));

        max_error = max([max(max(abs(V_c) - V_max), 0), max(max(V_min - abs(V_c)),0), max_P_test, max_Q_test]);
        if max_error < error_save
            max_u = max(max(abs(V_c) - V_max), 0);
            max_l = max(max(V_min - abs(V_c)),0);
            max_P = max_P_test;
            max_Q = max_Q_test;
        end

        V = V_c;
        curtailment_cost = sum((max(real(V .* (conj(Y) * conj(V))) - P_D,0) + max(imag(V .* (conj(Y) * conj(V))) - Q_D,0))) * param.curtailment_cost_plus - sum((min(real(V .* (conj(Y) * conj(V))) - P_D,0) + min(imag(V .* (conj(Y) * conj(V))) - Q_D,0))) * param.curtailment_cost_minus;
        fprintf(param.log_file, 'Correction with fminunc. Inital error: %g, final error: %g, curtailment cost: %f\n', error_save, max_error, curtailment_cost); 
    end
end

function error = mean_deviation(V, P_D, Q_D, Y)
    V = complex(V(1:length(V)/2), V(length(V)/2+1:end));

    P_V = real(V .* (conj(Y) * conj(V)));
    Q_V = imag(V .* (conj(Y) * conj(V)));
    error = norm(P_V-P_D, 1.25)^1.25 + norm(Q_V-Q_D, 1.25)^1.25;
end

function [c, ceq] = voltage_constraints(V, V_max, V_min)
    V = complex(V(1:length(V)/2), V(length(V)/2+1:end));
    max_u = norm( max(abs(V) - (V_max), 0), 2)^2;
    max_l = norm( max((V_min) - abs(V), 0), 2)^2;
    c = max_u + max_l;
    ceq = [];
end
