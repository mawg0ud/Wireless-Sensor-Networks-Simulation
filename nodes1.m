function [pij, node_degree] = nodes1(cen1, neighbor_limit, alpha, Band, lamdamV, Thr)
    % Calculate distance matrix
    dij = dist(cen1');
    
    % Initialize power matrix
    pij = zeros(size(cen1, 1));

    % Initialize node degree vector
    node_degree = zeros(size(cen1, 1), 1);

    % Loop over each station
    for k = 1:size(cen1, 1)
        % Calculate neighbors within range for each station
        neighbors = find((dij(k, :) > 0) & (dij(k, :) <= neighbor_limit));

        % Update node degree
        node_degree(k) = length(neighbors);

        % Calculate power for each neighbor
        for neighbor_idx = 1:length(neighbors)
            neighbor = neighbors(neighbor_idx);
            theta = theta1(cen1(k, :), cen1(neighbor, :));
            pij(k, neighbor) = (dij(k, neighbor) ^ alpha) * theta;
        end
    end
end
