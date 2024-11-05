close all
clear
clc

%% Figure 4.1 - Position
phi = linspace(0, 2*pi, 10);

theta = [0, pi/4, pi/2, 3*pi/4, pi];

[PHI, THETA] = ndgrid(phi, theta);

PHI = reshape(PHI, [], 1);
THETA = reshape(THETA, [], 1);

ball = [PHI, THETA].';

ball = [ball, [0;0],[0;pi]];

points = [];

r = 0.05;
for i = 1:size(ball, 2)
    phi = ball(1, i);
    theta = ball(2, i);
    x = r*sin(theta)*cos(phi);
    y = r*sin(theta)*sin(phi);
    z = r*cos(theta);

    points = [points, [x;y;z]];
end

x = points(1,:);
y = points(2,:);
z = points(3,:) + 0.15;

figure;
hold on

plot3(x, y, z, 'o'); % 'o' specifies circle markers

grid on; % Adds grid lines for better visualization
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
axis equal
title('3D Scatter Plot of Points');

% Create a sphere
[x, y, z] = sphere(50);  % 50 is the resolution of the sphere

% Scale the sphere to the desired radius and shift its center
radius = 0.05;
x = radius * x;
y = radius * y;
z = radius * z + 0.15;

% Plot the sphere
% figure
hold on
grid on
surf(x, y, z, 'FaceAlpha', 0.3, 'EdgeColor', 'none');  % 0.3 makes the sphere 30% opaque

% Set axis properties
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
% % title('Transparent Sphere with Center at (0, 0.1) and Radius 0.05');
% 
xlim([-0.075,0.075])
ylim([-0.075 0.075])
zlim([0 0.25])


%% Figure 4.2 - Orientation
figure
hold on
grid on;
axis equal

theta = [-pi/2, -pi/4, 0, pi/4, pi/2];
phi = [-pi/2, -pi/4, 0, pi/4, pi/2];

psi = [0];

[Theta, Phi, Psi] = ndgrid(theta, phi, psi);

Theta = reshape(Theta, [], 1);
Phi = reshape(Phi, [], 1);
Psi = reshape(Psi, [], 1);

Angles = [Theta, Phi, Psi].';


v = [0;0;1];

for i = 1:size(Angles, 2)
    theta = Angles(1,i);
    phi = Angles(2,i);
    psi = Angles(3,i);

    Rx_1 = [1 0 0;                  % rotate about x-axis
    0 cos(theta) -sin(theta);
    0 sin(theta) cos(theta)];
    
    Ry = [cos(phi) 0 sin(phi);    % rotate about y-axis
    0 1 0;
    -sin(phi) 0 cos(phi)];
    
    Rx_2 = [1 0 0;                  % rotate about x-axis
    0 cos(psi) -sin(psi);
    0 sin(psi) cos(psi)];

    R_star = Rx_1 * Ry * Rx_2;

    v_rotated = R_star * v *0.02;

    quiver3(0, 0, 0, v_rotated(1), v_rotated(2), v_rotated(3), 'LineWidth', 2, 'MaxHeadSize', 0.5, 'Color', 'r');
end
xlabel("X" ,'interpreter' , 'latex')
ylabel("Y" ,'interpreter' , 'latex')
zlabel("Z" ,'interpreter' , 'latex')