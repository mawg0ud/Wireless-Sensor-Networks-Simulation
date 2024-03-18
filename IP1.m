function [x, minpower, status] = IP1(CC, A_L_eq, B_L_eq, A_eq, B_eq, lb, ub, M, e)
% IP1: Solves the mixed integer linear programming problem.
% Inputs:
%   - CC: Cost coefficients for the objective function.
%   - A_L_eq, B_L_eq: Matrices defining the less than or equal constraints.
%   - A_eq, B_eq: Matrices defining the equality constraints.
%   - lb, ub: Lower and upper bounds for decision variables.
%   - M: Set of indices for binary variables.
%   - e: Small value to avoid numerical issues.
% Outputs:
%   - x: Solution vector.
%   - minpower: Minimum power achieved.
%   - status: Status of the optimization (1 if successful, 0 otherwise).

% Define objective function
f = CC';

% Set optimization options
opts = optimoptions('intlinprog', 'Display', 'off');

% Solve the mixed integer linear programming problem
[x, minpower, status] = intlinprog(f, M, A_L_eq, B_L_eq, A_eq, B_eq, lb, ub, opts);

% Check if the solution is valid
if status ~= 1
    x = [];
    minpower = NaN;
end

end
