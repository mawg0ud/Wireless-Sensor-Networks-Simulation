
% Load data from fixedTopology_n20_m30_lim15.mat file
load fixedTopology_n20_m30_lim15

% Set the neighbor limit
neighbor_limit = 25;

% Calculate distance matrix
dij = dist(cen1');

% Initialize variables
alpha = 2;
Band = 2000;
lamdamV = 0.06 * Band:0.06 * Band:0.9 * Band;
Thr = 0:100:300;

% Loop over different Thr values
for vv = 1:3
    % Loop over different lamdam values
    for lamdam = 0.02 * Band:0.02 * Band:0.9 * Band
        % Perform computations for each lamdam value
        disp(['Current lamdam: ', num2str(lamdam)]);
    end
end
