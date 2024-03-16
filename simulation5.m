

clear all
close all
clc
warning off  
echo off
format short g

% Assume minimum distance in X, and Y between stations
string_table = [];

% Assume number of stations n
n = 20;

% All the nodes with zero power
% pick any one to consider it max
imax = randi(n);

load fixedTopology_n20_m30_lim15
neighbor_limit = 25;
dij = dist(cen1');

% assuming a value to alpha
alpha = 2;

Band = 2000;

lamdamV = 0.06 * Band:0.06 * Band;

% Threshold values
Thr = 0:100:300;
Thr_count = 0;

for vv = 1:3
    Thr_count = Thr_count + 1;
    countLamda = 0;
    counter = 0;

    for lamdam = 0.02 * Band:0.02 * Band
        countLamda = countLamda + 1;

        % the power matrix
        pij = zeros(n);

        for n_runs = 1:10
            % Generate random numbers with mean lamdam and variance lamdam/2
            % Pick a request lamda_s_d

            % Calculate stations within range to each one
            [pij, node_degree] = nodes1(cen1, neighbor_limit, alpha, Band, lamdam, Thr);

            % Fixed Number of requests for all runs for threshold analysis
            load Req_Bank_n_20_B_500_small_lamdam;
            reqmatAll = eval(['reqmat_run_', num2str(n_runs)]);
            Ind1 = find(reqmatAll(:, 4) == lamdam);
            reqmat = reqmatAll(Ind1, 1:3);

            % Initilization
            consP = 0;
            Tlamda_sd = 0;

            % consider number of requests
            for nr = 1:size(reqmat, 1)
                lamda_sd = reqmat(nr, 1);
                s = reqmat(nr, 2);
                d = reqmat(nr, 3);

                % Generate the power vector
                PVij = transpose(dij .^ alpha);
                Vxij = ones(n);  % Assuming all nodes are within range
                PVij = PVij(:)';
                Vxij = Vxij(:)';

                CC = PVij;

                % Power constraint
                PPmax = 3000;
                ppi = sum(pij');
                max_node = find(ppi == max(ppi));
                max_node = max_node(1);

                if max(ppi) > 0
                    TThr = Thr(vv);
                    P_constraint = mean(ppi(find(ppi))) + Thr(vv);
                else
                    P_constraint = Thr(vv);
                end

                A_L_eq = [];
                B_L_eq = [];

                % Power constraint as in the old paper
                for i = 1:n
                    pci = [zeros(length(1:(i-1)), n); pij(i, :) + dij(i, :) .^ alpha; zeros(length((i+1):n), n)];
                    pci = pci';
                    pci = pci(:);
                    cond1 = (ppi(s) < P_constraint) && (max(ppi) > 0);
                    if cond1
                        A_L_eq = [A_L_eq; pci];
                        B_L_eq = [B_L_eq; P_constraint];
                    end

                    for j = 1:n
                        A_L_eq = [A_L_eq; pci];
                        ddij = PPmax - (dij(i, j) .^ alpha) * Vxij(i, j);
                        B_L_eq = [B_L_eq; ddij];
                        A_L_eq = [A_L_eq; pci];
                        B_L_eq = [B_L_eq; PPmax];
                    end
                end

                % Delay constraint
                delta_sd = ceil(2 * n / 3);
                A_L_eq = [A_L_eq; ones(1, n*n)];
                B_L_eq = [B_L_eq; delta_sd];

                % Bandwidth constraint
                AAA = [];
                for i = 1:n
                    AA = zeros(n);
                    AA(i, :) = 1;
                    AA(:, i) = AA(:, i) + 1;
                    AA = AA';
                    AA = AA(:)';
                    AAA = [AAA; AA];
                end

                A_L_eq = [A_L_eq; AAA];
                B_L_eq = [B_L_eq; (Band/lamda_sd) * ones(n, 1)];

                % Solve optimization problem
                lb = zeros(1, n*n);
                ub = ones(1, n*n);
                e = 2^-24;
                M = 1:n*n;

                [x, minpower, status] = IP1(CC, A_L_eq, B_L_eq, lb, ub, M, e);

                if status == 1
                    counter = counter + 1;
                    consP = consP + minpower;
                    pij = pij + reshape(x, n, n) .* (dij .^ alpha);
                    Tlamda_sd = Tlamda_sd + lamda_sd;
                end
            end

            % Store results
            ppi = sum(pij');
            eval(['status_Thr_', num2str(Thr_count), '_run_', num2str(n_runs), ...
                  '_lamda_', num2str(countLamda), '= status;']);
            if status == 1
                % Update string table
                % Update other metrics as needed
            end
        end
    end
end
