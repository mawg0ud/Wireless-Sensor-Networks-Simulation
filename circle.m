function [x, y] = circle(center, radius, num_points)
% circle: Generates points on a circle.
% Inputs:
%   - center: Coordinates of the center of the circle [x, y].
%   - radius: Radius of the circle.
%   - num_points: Number of points to generate (optional, default is 100).
% Outputs:
%   - x: x-coordinates of the points on the circle.
%   - y: y-coordinates of the points on the circle.

if nargin < 3
    num_points = 100; % Default number of points
end

theta = linspace(0, 2*pi, num_points);
x = center(1) + radius * cos(theta);
y = center(2) + radius * sin(theta);

end
