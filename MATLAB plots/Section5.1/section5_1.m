close all
clear
clc

%% Figure 5.1 - Orientation
figure
hold on
grid on
axis equal
theta = linspace(0, 2*pi, 9);
theta(end) = [];
phi = linspace(-pi/2, pi/2, 5);
phi(end) = [];
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

    quiver3(0, 0, 0, v_rotated(1), v_rotated(2), v_rotated(3), 'LineWidth', 2, 'MaxHeadSize', 0.5, 'Color', 'b');
end
xlabel("X" ,'interpreter' , 'latex')
ylabel("Y" ,'interpreter' , 'latex')
zlabel("Z" ,'interpreter' , 'latex')