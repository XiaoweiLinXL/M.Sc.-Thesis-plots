close all
clear
clc

%% Figure 1.3 - Workspace
% Create a sphere
[x, y, z] = sphere(50);  % 50 is the resolution of the sphere

% Scale the sphere to the desired radius and shift its center
radius = 0.05;
x = radius * x;
y = radius * y;
z = radius * z + 0.15;

% Plot the sphere
figure
hold on
grid on
surf(x, y, z, 'FaceAlpha', 0.3, 'EdgeColor', 'none');  % 0.3 makes the sphere 30% opaque

% Set axis properties
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');

xlim([-0.075,0.075])
ylim([-0.075 0.075])
zlim([0 0.25])


drawCylinder([0,0.05,0], [1;0;0;0], 0.5e-2, 1e-2)
drawCylinder([0.05,0,0], [1;0;0;0], 0.5e-2, 1e-2)
drawCylinder([0,-0.05,0], [1;0;0;0], 0.5e-2, 1e-2)
drawCylinder([-0.05,0,0], [1;0;0;0], 0.5e-2, 1e-2)