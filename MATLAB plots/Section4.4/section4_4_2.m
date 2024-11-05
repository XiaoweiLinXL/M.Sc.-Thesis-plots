close all
clear
clc

%%
% Spheres
unit = 'decimeter';
if strcmp(unit, 'meter') 
    scale = 1;
elseif strcmp(unit, 'decimeter')
    scale = 0.1;
elseif strcmp(unit, 'centimeter')
    scale = 0.01;
elseif strcmp(unit, 'millimeter')
    scale = 0.001;
end

mu0 =  4*pi*1e-7; % air permeability
cylinder_dia = 6.35e-3 / scale; % cylinder diameter 6.35e-3 m (1/4 inch)
                                % converted to corresponding unit
                                % using the scale factor
cylinder_length = 12.7e-3 /scale; % cylinder length 12.7e-3 m

Volume = cylinder_length*pi*(norm(cylinder_dia)/2)^2; % Volume of the sphere dipole
                                                      % in the unit of (specified unit)^3

B_r = 1.32; % Residual flux density (T)
                            
sph_dip = B_r*Volume/mu0; % spheres magnetic dipole mu_norm = B_r*V/mu0

type = "3D"; % 3-axis sensor

%% Magnet configurations
% Position using spherical coordinate
phi_sphere = linspace(0, 2*pi, 10);
theta_sphere = [pi/4, pi/2, 3*pi/4];

% Orientation
theta = [-pi/2, -pi/4, 0, pi/4, pi/2];
phi = [-pi/2, -pi/4, 0, pi/4, pi/2];
psi = [0];

[PHI_sphere, THETA_sphere, Theta, Phi, Psi] = ndgrid(phi_sphere, theta_sphere, theta, phi, psi);

PHI_sphere = reshape(PHI_sphere, [], 1);
THETA_sphere = reshape(THETA_sphere, [], 1);
Theta = reshape(Theta, [], 1);
Phi = reshape(Phi, [], 1);
Psi = reshape(Psi, [], 1);

sphere = [PHI_sphere, THETA_sphere, Theta, Phi, Psi].';

phi_sphere_extreme = [0];
theta_sphere_extreme = [0,pi];
[PHI_sphere, THETA_sphere, Theta, Phi, Psi] = ndgrid(phi_sphere_extreme, theta_sphere_extreme, theta, phi, psi);
PHI_sphere = reshape(PHI_sphere, [], 1);
THETA_sphere = reshape(THETA_sphere, [], 1);
Theta = reshape(Theta, [], 1);
Phi = reshape(Phi, [], 1);
Psi = reshape(Psi, [], 1);

sphere = [sphere, [PHI_sphere, THETA_sphere, Theta, Phi, Psi].'];

points = [];

r = 5*0.01; % 5*sqrt(3) cm
for i = 1:size(sphere, 2)
    phi = sphere(1, i);
    theta = sphere(2, i);
    x = r*sin(theta)*cos(phi)/scale;
    y = r*sin(theta)*sin(phi)/scale;
    z = (r*cos(theta)+0.15)/scale;

    points = [points, [x;y;z;sphere(3,i);sphere(4,i);sphere(5,i)]];
end

magnet_conf = points;

%% 3 sensor on a circle
sens_num = 3;
[rcond_all_sens_conf_circle, mean_min_svd_all_sens_conf_circle] = ...
    criteria_from_sensor_config(sens_num, scale, magnet_conf, B_r, Volume, type);

% Pareto front 
figure
plot_points_with_arrows(rcond_all_sens_conf_circle, mean_min_svd_all_sens_conf_circle, [0.8500, 0.3250, 0.0980], 0.2, 1.2e-5)
hold on
 
grid on

load('results_3_3_axis_multiobj_1000gen_1000pop_workspace048_distance10_finer_smaller_sphere_bigger_workspace_mean_sigma_min.mat')
scatter(-fval(:,1),-fval(:,2), 10, 'MarkerEdgeColor', [0.4660, 0.6740, 0.1880], 'MarkerFaceColor',[0.4660, 0.6740, 0.1880]);
hold on
xlabel('$\min_{x\in \mathcal{X}} \frac{1}{\kappa}(x)$', 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
ylabel('$mean_{x \in \mathcal{X}} \sigma_{min}(x)$', 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
xlim([0 0.2])
ylim([0 1.2e-5])

legend('sensors on circle', 'pareto front form GA', 'Location', 'southeast')

%% 3 sensor on a circle - visualize configuration
load('results_3_3_axis_multiobj_1000gen_1000pop_workspace048_distance10_finer_smaller_sphere_bigger_workspace_mean_sigma_min.mat')

sigma_min = -fval(:,2);

sorted_sigma_min = sort(sigma_min);
start_point_sigma_min = sorted_sigma_min(20);
mid_point_sigma_min = sorted_sigma_min(length(sorted_sigma_min)/2);
end_point_sigma_min = sorted_sigma_min(end);

index = [find(abs(sigma_min-start_point_sigma_min)<1e-10), 
         find(abs(sigma_min-mid_point_sigma_min)<1e-10),
         find(abs(sigma_min-end_point_sigma_min)<1e-10)];

for i = 1:length(index)
    sol_selected = sol(index(i),:);
    
    % Plot sensor config
    plot_sensor_config(sol_selected);
end

%% 4 sensor on a circle
sens_num = 4;
[rcond_all_sens_conf_circle, mean_min_svd_all_sens_conf_circle] = ...
    criteria_from_sensor_config(sens_num, scale, magnet_conf, B_r, Volume, type);

% Pareto front 
figure
plot_points_with_arrows(rcond_all_sens_conf_circle, mean_min_svd_all_sens_conf_circle, [0.8500, 0.3250, 0.0980], 0.25, 1.5e-5)
hold on
 
grid on

load('results_4_3_axis_multiobj_1000gen_1000pop_workspace048_distance10_finer_smaller_sphere_bigger_workspace_mean_sigma_min.mat')
scatter(-fval(:,1),-fval(:,2), 10, 'MarkerEdgeColor', [0.4660, 0.6740, 0.1880], 'MarkerFaceColor',[0.4660, 0.6740, 0.1880]);
hold on
xlabel('$\min_{x\in \mathcal{X}} \frac{1}{\kappa}(x)$', 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
ylabel('$mean_{x \in \mathcal{X}} \sigma_{min}(x)$', 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
xlim([0 0.25])
ylim([0 1.5e-5])

legend('sensors on circle', 'pareto front form GA', 'Location', 'southeast')

%% 4 sensor on a circle - visualize configuration
load('results_4_3_axis_multiobj_1000gen_1000pop_workspace048_distance10_finer_smaller_sphere_bigger_workspace_mean_sigma_min.mat')

sigma_min = -fval(:,2);

sorted_sigma_min = sort(sigma_min);
start_point_sigma_min = sorted_sigma_min(20);
mid_point_sigma_min = sorted_sigma_min(length(sorted_sigma_min)/2);
end_point_sigma_min = sorted_sigma_min(end);

index = [find(abs(sigma_min-start_point_sigma_min)<1e-10), 
         find(abs(sigma_min-mid_point_sigma_min)<1e-10),
         find(abs(sigma_min-end_point_sigma_min)<1e-10)];

for i = 1:length(index)
    sol_selected = sol(index(i),:);
    
    % Plot sensor config
    plot_sensor_config(sol_selected);
end

%% 5 sensor on a circle
sens_num = 5;
[rcond_all_sens_conf_circle, mean_min_svd_all_sens_conf_circle] = ...
    criteria_from_sensor_config(sens_num, scale, magnet_conf, B_r, Volume, type);

% Pareto front 
figure
plot_points_with_arrows(rcond_all_sens_conf_circle, mean_min_svd_all_sens_conf_circle, [0.8500, 0.3250, 0.0980], 0.3, 1.8e-5)
hold on
 
grid on

load('results_5_3_axis_multiobj_1000gen_1000pop_workspace048_distance10_finer_smaller_sphere_bigger_workspace_mean_sigma_min.mat')
scatter(-fval(:,1),-fval(:,2), 10, 'MarkerEdgeColor', [0.4660, 0.6740, 0.1880], 'MarkerFaceColor',[0.4660, 0.6740, 0.1880]);
hold on
xlabel('$\min_{x\in \mathcal{X}} \frac{1}{\kappa}(x)$', 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
ylabel('$mean_{x \in \mathcal{X}} \sigma_{min}(x)$', 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
xlim([0 0.30])
ylim([0 1.8e-5])

legend('sensors on circle', 'pareto front form GA', 'Location', 'southeast')

%% 5 sensor on a circle - visualize configuration
load('results_5_3_axis_multiobj_1000gen_1000pop_workspace048_distance10_finer_smaller_sphere_bigger_workspace_mean_sigma_min.mat')

sigma_min = -fval(:,2);

sorted_sigma_min = sort(sigma_min);
start_point_sigma_min = sorted_sigma_min(20);
mid_point_sigma_min = sorted_sigma_min(length(sorted_sigma_min)/2);
end_point_sigma_min = sorted_sigma_min(end);

index = [find(abs(sigma_min-start_point_sigma_min)<1e-10), 
         find(abs(sigma_min-mid_point_sigma_min)<1e-10),
         find(abs(sigma_min-end_point_sigma_min)<1e-10)];

for i = 1:length(index)
    sol_selected = sol(index(i),:);
    
    % Plot sensor config
    plot_sensor_config(sol_selected);
end

%% 9 sensor on a circle
sens_num = 9;
[rcond_all_sens_conf_circle, mean_min_svd_all_sens_conf_circle] = ...
    criteria_from_sensor_config(sens_num, scale, magnet_conf, B_r, Volume, type);

% Pareto front 
figure
plot_points_with_arrows(rcond_all_sens_conf_circle, mean_min_svd_all_sens_conf_circle, [0.8500, 0.3250, 0.0980], 0.3, 2.5e-5)
hold on
 
grid on

load('results_9_3_axis_multiobj_1000gen_1000pop_workspace048_distance10_finer_smaller_sphere_bigger_workspace_mean_sigma_min.mat')
scatter(-fval(:,1),-fval(:,2), 10, 'MarkerEdgeColor', [0.4660, 0.6740, 0.1880], 'MarkerFaceColor',[0.4660, 0.6740, 0.1880]);
hold on
xlabel('$\min_{x\in \mathcal{X}} \frac{1}{\kappa}(x)$', 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
ylabel('$mean_{x \in \mathcal{X}} \sigma_{min}(x)$', 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
xlim([0 0.30])
ylim([0 2.5e-5])

legend('sensors on circle', 'pareto front form GA', 'Location', 'southeast')

%% 9 sensor on a circle - visualize configuration
load('results_9_3_axis_multiobj_1000gen_1000pop_workspace048_distance10_finer_smaller_sphere_bigger_workspace_mean_sigma_min.mat')

sigma_min = -fval(:,2);

sorted_sigma_min = sort(sigma_min);
start_point_sigma_min = sorted_sigma_min(10);
mid_point_sigma_min = sorted_sigma_min(length(sorted_sigma_min)/2);
end_point_sigma_min = sorted_sigma_min(end);

index = [find(abs(sigma_min-start_point_sigma_min)<1e-10), 
         find(abs(sigma_min-mid_point_sigma_min)<1e-10),
         find(abs(sigma_min-end_point_sigma_min)<1e-10)];

for i = 1:length(index)
    sol_selected = sol(index(i),:);
    
    % Plot sensor config
    plot_sensor_config(sol_selected);
end

%% 16 sensor on a circle
sens_num = 16;
[rcond_all_sens_conf_circle, mean_min_svd_all_sens_conf_circle] = ...
    criteria_from_sensor_config(sens_num, scale, magnet_conf, B_r, Volume, type);

% Pareto front 
figure
plot_points_with_arrows(rcond_all_sens_conf_circle, mean_min_svd_all_sens_conf_circle, [0.8500, 0.3250, 0.0980], 0.3, 3.0e-5)
hold on
 
grid on

load('results_16_3_axis_multiobj_1000gen_1000pop_workspace048_distance10_finer_smaller_sphere_bigger_workspace_mean_sigma_min.mat')
scatter(-fval(:,1),-fval(:,2), 10, 'MarkerEdgeColor', [0.4660, 0.6740, 0.1880], 'MarkerFaceColor',[0.4660, 0.6740, 0.1880]);
hold on
xlabel('$\min_{x\in \mathcal{X}} \frac{1}{\kappa}(x)$', 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
ylabel('$mean_{x \in \mathcal{X}} \sigma_{min}(x)$', 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
xlim([0 0.30])
ylim([0 3.0e-5])

legend('sensors on circle', 'pareto front form GA', 'Location', 'southeast')

%% 16 sensor on a circle - visualize configuration
load('results_16_3_axis_multiobj_1000gen_1000pop_workspace048_distance10_finer_smaller_sphere_bigger_workspace_mean_sigma_min.mat')

sigma_min = -fval(:,2);

sorted_sigma_min = sort(sigma_min);
start_point_sigma_min = sorted_sigma_min(10);
mid_point_sigma_min = sorted_sigma_min(length(sorted_sigma_min)/2);
end_point_sigma_min = sorted_sigma_min(end);

index = [find(abs(sigma_min-start_point_sigma_min)<1e-10), 
         find(abs(sigma_min-mid_point_sigma_min)<1e-10),
         find(abs(sigma_min-end_point_sigma_min)<1e-10)];

for i = 1:length(index)
    sol_selected = sol(index(i),:);
    
    % Plot sensor config
    plot_sensor_config(sol_selected);
end
%% Sensor number vs best achievable criteria
% Figure 4.20 - 4.23
sensor_number = 3:1:20;

best_mean_sigma_min = [1.10e-5, 1.38e-5, 1.56e-5, 1.72e-5, 1.85e-5, 1.98e-5, 2.10e-5, 2.21e-5, 2.32e-5, 2.42e-5, ...
                       2.52e-5, 2.62e-5, 2.71e-5, 2.80e-5, 2.89e-5, 2.97e-5, 3.05e-5, 3.13e-5];
best_mean_sigma_min_diameter = [0.069,0.072,0.076,0.076,0.076,0.076,0.076,0.076,0.076,0.076,0.076,0.076,0.076,0.076,0.076,0.076,0.076,0.076]*2;

best_min_rcond = [0.17175,0.243515,0.254841,0.254841,0.254841,0.254841,0.254841,0.254841,0.254841,0.254841,...
                  0.254841,0.254841,0.254841,0.254841,0.254841,0.254841,0.254841,0.254841];
best_min_rcond_diameter = [0.378,0.355,0.329,0.329,0.329,0.329,0.329,0.329,0.329,0.329,0.329,0.329,0.329,0.329,0.329,0.329,0.329,0.329]*2;


figure
plot(sensor_number, best_mean_sigma_min, LineWidth=2);
xlabel("Sensor number", 'Interpreter', 'latex');
ylabel("Best $mean_{x \in \mathcal{X}} \sigma_{min}(x)$", 'Interpreter', 'latex');
grid on

figure
plot(sensor_number, best_mean_sigma_min_diameter, LineWidth=2);
xlabel("Sensor number", 'Interpreter', 'latex');
ylabel("Circle diameter (m)", 'Interpreter', 'latex');
grid on


figure 
plot(sensor_number, best_min_rcond, LineWidth=2);
xlabel("Sensor number", 'Interpreter', 'latex');
ylabel("Best $\min_{x \in \mathcal{X} \frac{1}{\kappa}(x)}$", 'Interpreter', 'latex');
grid on

figure
plot(sensor_number, best_min_rcond_diameter, LineWidth=2);
xlabel("Sensor number", 'Interpreter', 'latex');
ylabel("Circle diameter (m)", 'Interpreter', 'latex');
grid on


%% sigma_min to 998 percentile noise
% Figure 4.24 - 4.25
number_of_sensor_range = 3:1:20;
percentile998_noise_norm = [];

% Chi distribution
for index = 1:length(number_of_sensor_range)
    number_of_sensor = number_of_sensor_range(index);
    axes = 3*number_of_sensor;
    noise_std_one_axis = 5e-08;
    noise_variance_one_axis = noise_std_one_axis^2;
    noise_variance_measurement_pair = 2*noise_variance_one_axis;
    
    x = 0:1e-10:5e-6;
    
    cdf_analytical = [];
    for i = 1:length(x)
        cdf_analytical_onepoint = gammainc(axes/2,(x(i)^2)/(2*noise_variance_measurement_pair));
        cdf_analytical = [cdf_analytical, cdf_analytical_onepoint];
    end
    cdf_analytical=1-cdf_analytical;
        
    for j = 1:length(x)
        if cdf_analytical(j) >= 0.99
            threshold = x(j);
            break
        end
    end

    percentile998_noise_norm = [percentile998_noise_norm, threshold];
end

sensor_number = 3:1:20;

% min(sigma_min) when the sensor configuration has the best mean(sigma_min)
min_sigma_min = [1.42e-6, 1.68e-6, 1.92e-6, 2.11e-6, 2.28e-6, 2.43e-6, 2.58e-6, 2.72e-6, 2.85e-6, 2.98e-6, ...
                         3.10e-6, 3.22e-6, 3.33e-6, 3.44e-6, 3.55e-6, 3.65e-6, 3.75e-6, 3.85e-6];

signal_to_noise_ratio = min_sigma_min./percentile998_noise_norm;


figure
hold on
plot(sensor_number, min_sigma_min, LineWidth=2);
plot(sensor_number, percentile998_noise_norm, LineWidth=2);
xlabel("Sensor number", 'Interpreter', 'latex');
ylabel("Signal strength (T)", 'Interpreter', 'latex');
grid on

figure
plot(sensor_number, signal_to_noise_ratio, LineWidth=2);
xlabel("Sensor number", 'Interpreter', 'latex');
ylabel("Signal to noise ratio", 'Interpreter', 'latex');
grid on

%% General function
% Quaternion to rotation matrix
function R = quaternionToMatrix(q)
    % Extract the scalar and vector parts from the quaternion
    w = q(1);
    x = q(2);
    y = q(3);
    z = q(4);

    % Calculate the elements of the rotation matrix
    R11 = 1 - 2*y^2 - 2*z^2;
    R12 = 2*x*y - 2*z*w;
    R13 = 2*x*z + 2*y*w;
    R21 = 2*x*y + 2*z*w;
    R22 = 1 - 2*x^2 - 2*z^2;
    R23 = 2*y*z - 2*x*w;
    R31 = 2*x*z - 2*y*w;
    R32 = 2*y*z + 2*x*w;
    R33 = 1 - 2*x^2 - 2*y^2;

    % Combine the elements into the rotation matrix
    R = [R11, R12, R13;
         R21, R22, R23;
         R31, R32, R33];
end

% Convert a vector to a skew symmetric matrix
function v_hat = skew(v)
    v_hat = [0 -v(3) v(2); 
             v(3) 0 -v(1);
             -v(2) v(1) 0];
end

%% Criteria from sensor config
function [rcond_all_sens_conf_circle, mean_min_svd_all_sens_conf_circle] = ...
    criteria_from_sensor_config(sens_num, scale, magnet_conf, B_r, Volume, type)
    % Define the radius
    radius = 0.01:0.001:1;
    
    % Initialize a cell array to hold all coordinates for each value of k
    all_coordinates = cell(1, length(radius));
    
    for i = 1:length(radius)
        % Angles for the points on the circle
        angles = linspace(0, 2*pi, sens_num+1);
        angles(end) = []; % Remove the last element to get k points
        
        % Calculate coordinates for the circle
        x = radius(i) * cos(angles);
        y = radius(i) * sin(angles);
        
        % Store the coordinates
        all_coordinates{i} = [x', y'].';
    end
    
    % Iterate throught the sensor configs to get the pareto front
    rcond_all_sens_conf_circle = [];
    mean_min_svd_all_sens_conf_circle = [];
    
    for i = 1:length(all_coordinates)
        sens_pos = all_coordinates{i};
        sens_num = size(sens_pos,2);
    
        % Add the line of zero as the z coordinate to the position of sensors
        z_coordinate = zeros(1, size(sens_pos, 2));
        sens_pos = [sens_pos; z_coordinate];
        sens_pos = sens_pos/scale;
    
        % Default orientation
        default_or = [1;0;0;0];
    
        % Orientation for all sensors
        sens_or_unitary = repmat(default_or, 1, sens_num);
        
        rcond_one_sens_conf = [];
        min_svd_one_sens_conf = [];
        average_sigma_one_sens_conf = [];
    
        % Collect the reciprocal condition number for each magnet configuration
        for magnet_num=1:size(magnet_conf,2)
            magnet_pos = magnet_conf(1:3,magnet_num);
            theta = magnet_conf(4, magnet_num);
            phi = magnet_conf(5, magnet_num);
            psi = magnet_conf(6, magnet_num);
        
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
            R_star = Ry * Rx_1;
        
            J_scaled = [];
            for j=1:sens_num
                J_scaled = [J_scaled;J_analytical_sensor_B_to_world_pR(sens_pos(:,j), sens_or_unitary(:,j), ...
                    magnet_pos, B_r, Volume, R_star, type)];
            end
        
            sigma_scaled = svd(J_scaled);
        
            num_dof = 5;
            reciprocal_condition_number_scaled = sigma_scaled(num_dof)/sigma_scaled(1);
            
            rcond_one_sens_conf = [rcond_one_sens_conf; reciprocal_condition_number_scaled];
            min_svd_one_sens_conf = [min_svd_one_sens_conf; sigma_scaled(num_dof)];
    
            % Exclude the 0 sigma_min
            sigma_scaled(end) = [];
            average_sigma_one_sens_conf = [average_sigma_one_sens_conf; mean(sigma_scaled)];
        end
    
        index_min_rcond = find(rcond_one_sens_conf == min(rcond_one_sens_conf));
        index_min_minsvd = find(min_svd_one_sens_conf == min(min_svd_one_sens_conf));
    
        % Put the min value in the workspace into the list
        rcond_all_sens_conf_circle = [rcond_all_sens_conf_circle; min(rcond_one_sens_conf)];
        mean_min_svd_all_sens_conf_circle = [mean_min_svd_all_sens_conf_circle; mean(min_svd_one_sens_conf)];
    end
end

%% Plot sensor config
function [] = plot_sensor_config(sol)
    % Constants
    sens_dia = 2e-2;
    sens_hi = 1e-2;
    
    % Define the radius and center
    r = 0.05;
    centerX = 0;
    centerY = 0;
    
    % Parametric equation for the circle
    theta = linspace(0, 2*pi, 100); % Divide the circle into 100 points
    x = centerX + r * cos(theta);
    y = centerY + r * sin(theta);
    
    % Plot the circle
    figure
    plot(x, y, 'LineWidth', 2)
    grid on 
    hold on
    
    % Plot sensor config
    sens_conf = sol;
    sens_num = sens_conf(end);
    sens_conf(end) = [];
    
    sens_conf = sens_conf*0.1;
    sens_conf = reshape(sens_conf, 2, []);
    
    % Add the line of zero as the z coordinate to the position of sensors
    z_coordinate = zeros(1, size(sens_conf, 2));
    sens_pos = [sens_conf; z_coordinate];
    
    % Default orientation
    default_or = [1;0;0;0];
    
    % Orientation for all sensors
    sens_or = repmat(default_or, 1, sens_num);
    
    for i = 1:sens_num
        radius = sens_dia/2;
        centerX = sens_pos(1,i); % Replace with the x-coordinate of the center
        centerY = sens_pos(2,i); % Replace with the y-coordinate of the center
        
        % Parametric equation for the circle
        theta = linspace(0, 2*pi, 100); % Divide the circle into 100 points
        circleX = centerX + radius * cos(theta);
        circleY = centerY + radius * sin(theta);
        
        % Plot the circle
        plot(circleX, circleY, 'LineWidth', 2, 'Color', 'black')
    end
    
    view(0,90)
    legend("Magnet workspace" ,"Sensor")
    xlim([-0.5, 0.5])
    ylim([-0.5, 0.5])
    xlabel("X (m)", 'Interpreter', 'latex')
    ylabel("Y (m)", 'Interpreter', 'latex')
    axis equal
end



%% Jacobian analytical
% Jacobian analytical
% Jacobian of B in sensor frame to magnet position in world frame, p, and
% magnet orientation with respect to world frame, R.
function J = J_analytical_sensor_B_to_world_pR(sens_pos, sens_or, magnet_pos, B_r, Volume, R_star, type) % sens_pos, sens_or, magnet_pos, mu are all in world frame
    mu_world_hat = R_star * [0;0;1];

    % First part of the jacobian, from sensor reading to magnet position
    % Extract the rotation matrix of the sensor
    unit_sens_or = sens_or/norm(sens_or);
    sensor_rotation_matrix = quaternionToMatrix(unit_sens_or); % sensor to world frame

    % Calculate r in the world frame
    r_world = sens_pos - magnet_pos;
    r_world_hat = r_world/norm(r_world);
    
    % J of B with respect to p, evaluated at p*, R*
    J_position = -((3*B_r*Volume)/(4*pi*norm(r_world)^4)) * ((eye(3)-5*(r_world_hat*r_world_hat.')) * ...
                 (r_world_hat.'*mu_world_hat) + ...
                 mu_world_hat*r_world_hat.' + r_world_hat*mu_world_hat.');

    % Second part of the jacobian, from sensor reading to magnet
    % orientation R
    
    % J of B with respect to R, evaluated at p*, R*
    e3 = [0;0;1];
    J_angle = (B_r*Volume/(4*pi)) * (1/(norm(r_world)^3)) * (3*(r_world_hat*r_world_hat.')-eye(3)) * ...
              (-R_star*skew(e3));

    J = [J_position, J_angle];

    J = sensor_rotation_matrix.' * J;

    % If 1d sensor, select just the last row
    if type == "1D"
        J = J(3, :);
    end
end

%% Plot points with arrows
function [] = plot_points_with_arrows(x_coord, y_coord, color, xlim_max, ylim_max)
    side_length = 2 * 0.01:0.001:1;
    % Define the points
    points = [x_coord, y_coord];
    
    % Plot the points (for reference)
    plot(points(:,1), points(:,2), '-', 'Color',color, 'MarkerSize',2);
    hold on;
    
    % Set axis limits for better visualization
    xlim([0 xlim_max]);
    ylim([0 ylim_max]);
    
    % Get the current axes position
    ax = gca;
    ax.Units = 'normalized';
    ax_pos = ax.Position;
    
    % Normalize the coordinates to the figure
    x_range = ax.XLim(2) - ax.XLim(1);
    y_range = ax.YLim(2) - ax.YLim(1);
    
    % Normalize each point
    norm_points = zeros(size(points));
    for i = 1:size(points, 1)
        norm_points(i, 1) = (points(i, 1) - ax.XLim(1)) / x_range * ax_pos(3) + ax_pos(1);
        norm_points(i, 2) = (points(i, 2) - ax.YLim(1)) / y_range * ax_pos(4) + ax_pos(2);
    end
    
    % Add arrows at equal distances
    num_arrows = 10;
    step = 1;
    total_points = size(points, 1);
    interval = floor((total_points - 20) / (num_arrows + 1)); % Adjust interval calculation

    for i = 0:num_arrows-1
        start_idx = 20 + i * interval; % Start from the 10th coordinate
        end_idx = start_idx + step;
        if end_idx <= total_points
            arrow = annotation('arrow', [norm_points(start_idx, 1), norm_points(end_idx, 1)], ...
                                        [norm_points(start_idx, 2), norm_points(end_idx, 2)], ...
                                        'Color', color, 'LineWidth', 2); % Set the color and width of the arrow
        end
    end

    
    % Labels
    xlabel('X-axis');
    ylabel('Y-axis');
    hold off;
end