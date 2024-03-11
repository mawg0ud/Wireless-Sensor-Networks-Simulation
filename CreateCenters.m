

function [centers] = CreateCenters(num_centers, area_limits)
% CreateCenters: Generates random coordinates for sensor nodes within a specified area.
% Inputs:
%   - num_centers: Number of sensor nodes to generate.
%   - area_limits: [xmin, xmax, ymin, ymax] defining the bounding box for the sensor field.
% Output:
%   - centers: Matrix containing the coordinates of the generated sensor nodes.

% Generate random x and y coordinates within the specified area limits
x_coordinates = rand(1, num_centers) * (area_limits(2) - area_limits(1)) + area_limits(1);
y_coordinates = rand(1, num_centers) * (area_limits(4) - area_limits(3)) + area_limits(3);

% Combine x and y coordinates into a matrix
centers = [x_coordinates; y_coordinates]';

% Display the generated sensor nodes
figure;
scatter(centers(:,1), centers(:,2), 'filled');
xlabel('X-coordinate');
ylabel('Y-coordinate');
title('Generated Sensor Nodes');
axis([area_limits(1) area_limits(2) area_limits(3) area_limits(4)]);
axis equal;
grid on;

end
