close all
clear 
clc

%% Figure 3.3 - Positions
phi = linspace(0, 2*pi, 9);
phi(end) = [];
theta = linspace(0, pi, 7);
theta(1) = [];
theta(end) = [];
[PHI, THETA] = ndgrid(phi, theta);

PHI = reshape(PHI, [], 1);
THETA = reshape(THETA, [], 1);

ball = [PHI, THETA].';

% Top and bottom position
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


%% Figure 3.4 - Orientation cone facing up
figure
hold on
grid on;
axis equal

% Orientation
theta = linspace(-pi/2+deg2rad(30),pi/2-deg2rad(30),5);
phi = linspace(-pi/2+deg2rad(30),pi/2-deg2rad(30),5);
psi = 0;

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

    if i == (1+size(Angles, 2))/2
        quiver3(0, 0, 0, v_rotated(1), v_rotated(2), v_rotated(3), 'LineWidth', 2, 'MaxHeadSize', 0.5, 'Color', 'r');
    else
        quiver3(0, 0, 0, v_rotated(1), v_rotated(2), v_rotated(3), 'LineWidth', 2, 'MaxHeadSize', 0.5, 'Color', 'b');

    end

end
xlabel("X" ,'interpreter' , 'latex')
ylabel("Y" ,'interpreter' , 'latex')
zlabel("Z" ,'interpreter' , 'latex')

%% Figure 3.5 - Orientation cone facing down
figure
hold on
grid on;
axis equal

% Orientation
theta = linspace(pi/2+deg2rad(30),3*pi/2-deg2rad(30),5);
phi = linspace(-pi/2+deg2rad(30),pi/2-deg2rad(30),5);
psi = 0;

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

    if i == (1+size(Angles, 2))/2
        quiver3(0, 0, 0, v_rotated(1), v_rotated(2), v_rotated(3), 'LineWidth', 2, 'MaxHeadSize', 0.5, 'Color', 'r');
    else
        quiver3(0, 0, 0, v_rotated(1), v_rotated(2), v_rotated(3), 'LineWidth', 2, 'MaxHeadSize', 0.5, 'Color', 'b');

    end

end
xlabel("X" ,'interpreter' , 'latex')
ylabel("Y" ,'interpreter' , 'latex')
zlabel("Z" ,'interpreter' , 'latex')

%% Figure 3.5 - Orientation cone rotated 135 deg around y
figure
hold on
grid on;
axis equal

% Orientation
theta = linspace(-pi/2+deg2rad(30),pi/2-deg2rad(30),5);
phi = linspace(pi/4+deg2rad(30),5*pi/4-deg2rad(30),5);
psi = 0;

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

    % R_star = Rx_1 * Ry * Rx_2;
    R_star = Ry * Rx_1 * Rx_2;

    v_rotated = R_star * v *0.02;

    if i == (1+size(Angles, 2))/2
        quiver3(0, 0, 0, v_rotated(1), v_rotated(2), v_rotated(3), 'LineWidth', 2, 'MaxHeadSize', 0.5, 'Color', 'r');
    else
        quiver3(0, 0, 0, v_rotated(1), v_rotated(2), v_rotated(3), 'LineWidth', 2, 'MaxHeadSize', 0.5, 'Color', 'b');

    end

end
xlabel("X" ,'interpreter' , 'latex')
ylabel("Y" ,'interpreter' , 'latex')
zlabel("Z" ,'interpreter' , 'latex')