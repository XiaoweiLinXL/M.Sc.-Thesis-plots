%%
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


%% Grid sensor configuration
half_side_length = 0.16/scale;
LB = [-half_side_length, -half_side_length];
UB = [half_side_length, half_side_length];
z = 0;

number_of_sensor_on_side = 2:1:10;

sensor_grid_config = {};
for i=1:length(number_of_sensor_on_side)
    num_sample = number_of_sensor_on_side(i);
    
    C = cell(1, 2);
    [C{:}] = ndgrid(linspace(0, 1, num_sample));
    C = cellfun(@(a) a(:), C, 'Uni', 0);
    combo = [C{:}];
    conf = combo.*repmat(UB - LB, [size(combo, 1), 1]) + ...
        repmat(LB, [size(combo, 1), 1]);
    conf = conf.';
    
    conf = [conf;zeros(1,size(conf,2));ones(1,size(conf,2));zeros(1,size(conf,2));zeros(1,size(conf,2));zeros(1,size(conf,2))];

    sensor_grid_config{i} = conf;
end

% Compute the criteria
[reciprocal_condition_number_all_sensor_config, mean_min_singular_value_J_all_sensor_config] = ...
    criteria_from_sensor_config(sensor_grid_config, magnet_conf, B_r, Volume, type);

%% Figure 4.7, 4.8
number_of_sensor = [[2:1:10].^2];

figure
plot(number_of_sensor, reciprocal_condition_number_all_sensor_config, '-', 'Color', [0, 0.4470, 0.7410], "linewidth", 2)

xlabel('Number of sensors' , 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
ylabel('$\min_{x\in \mathcal{X}} \frac{1}{\kappa}(x)$', 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
ylim([0,0.3])
grid on

legend('grid 0.32m')

figure
plot(number_of_sensor, mean_min_singular_value_J_all_sensor_config, '-', 'Color', [0, 0.4470, 0.7410], "linewidth", 2)
legend('grid 0.32m')

xlabel('Number of sensors', 'FontSize', 12 , 'FontWeight', 'bold' ,'Interpreter', 'latex')
ylabel('$mean_{x \in \mathcal{X}} \sigma_{min}(x)$', 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
grid on


%% Grid sensor configuration with 5 sensor
half_side_length = 0.16/scale;
LB = [-half_side_length, -half_side_length];
UB = [half_side_length, half_side_length];
z = 0;

number_of_sensor_on_side = 2:1:10;

sensor_grid_config = {};
for i=1:length(number_of_sensor_on_side)
    num_sample = number_of_sensor_on_side(i);
    
    C = cell(1, 2);
    [C{:}] = ndgrid(linspace(0, 1, num_sample));
    C = cellfun(@(a) a(:), C, 'Uni', 0);
    combo = [C{:}];
    conf = combo.*repmat(UB - LB, [size(combo, 1), 1]) + ...
        repmat(LB, [size(combo, 1), 1]);
    conf = conf.';
    
    conf = [conf;zeros(1,size(conf,2));ones(1,size(conf,2));zeros(1,size(conf,2));zeros(1,size(conf,2));zeros(1,size(conf,2))];

    sensor_grid_config{i} = conf;
end

five_sensor_conf = [-half_side_length,-half_side_length,0,1,0,0,0;...
                    half_side_length,-half_side_length,0,1,0,0,0;...
                    -half_side_length,half_side_length,0,1,0,0,0;...
                    half_side_length,half_side_length,0,1,0,0,0;...
                    0,0,0,1,0,0,0].';

sensor_grid_config = [sensor_grid_config(1), {five_sensor_conf}, sensor_grid_config(2:end)];

% Compute the criteria
[reciprocal_condition_number_all_sensor_config, mean_min_singular_value_J_all_sensor_config] = ...
    criteria_from_sensor_config(sensor_grid_config, magnet_conf, B_r, Volume, type);

reciprocal_condition_number_all_sensor_config_1 = reciprocal_condition_number_all_sensor_config;
mean_min_singular_value_J_all_sensor_config_1 = mean_min_singular_value_J_all_sensor_config;

%% Figure 4.9
number_of_sensor = [[2:1:10].^2];
number_of_sensor = [number_of_sensor(1),5,number_of_sensor(2:end)];

figure
plot(number_of_sensor, reciprocal_condition_number_all_sensor_config_1, '-', 'Color', [0, 0.4470, 0.7410], "linewidth", 2)

xlabel('Number of sensors' , 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
ylabel('$\min_{x\in \mathcal{X}} \frac{1}{\kappa}(x)$', 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
ylim([0,0.3])
grid on

legend('grid 0.32m')


%% Square sensor configuration
half_side_length = 0.16/scale;
LB = [-half_side_length, -half_side_length];
UB = [half_side_length, half_side_length];
z = 0;

number_of_sensor_on_side = 2:1:26;

sensor_sqaure_config = {};
for i=1:length(number_of_sensor_on_side)
    num_sample = number_of_sensor_on_side(i);
    
    C = cell(1, 2);
    [C{:}] = ndgrid(linspace(0, 1, num_sample));
    C = cellfun(@(a) a(:), C, 'Uni', 0);
    combo = [C{:}];
    conf = combo.*repmat(UB - LB, [size(combo, 1), 1]) + ...
        repmat(LB, [size(combo, 1), 1]);
    conf = conf.';
    
    conf = [conf;zeros(1,size(conf,2));ones(1,size(conf,2));zeros(1,size(conf,2));zeros(1,size(conf,2));zeros(1,size(conf,2))];

    sensor_sqaure_config{i} = conf;
end

% Remove the sensors in the middle
for i = 1:length(sensor_sqaure_config)
    sensor_config = sensor_sqaure_config{i};

    % Initialize a logical array for columns to keep
    keep_columns = true(1, size(sensor_config, 2));

    % Loop through each column to check the conditions
    for j = 1:size(sensor_config, 2)
        if abs(sensor_config(1, j)) ~= half_side_length && abs(sensor_config(2, j)) ~= half_side_length
            keep_columns(j) = false;  % Mark column for deletion
        end
    end

    % Keep only the columns that meet the condition
    sensor_sqaure_config{i} = sensor_config(:, keep_columns);
end


% Compute the criteria
[reciprocal_condition_number_all_sensor_config, mean_min_singular_value_J_all_sensor_config] = ...
    criteria_from_sensor_config(sensor_sqaure_config, magnet_conf, B_r, Volume, type);

reciprocal_condition_number_all_sensor_config_2 = reciprocal_condition_number_all_sensor_config;
mean_min_singular_value_J_all_sensor_config_2 = mean_min_singular_value_J_all_sensor_config;



%% Sensors on circle
n = 4:100;

sensor_circle_config = {};
for i=1:length(n)
    theta = linspace(0, 2*pi, n(i)+1);
    theta(end) = [];

    x = half_side_length*sqrt(2) * cos(theta);
    y = half_side_length*sqrt(2) * sin(theta);
    conf = [x;y];
    conf = [conf;zeros(1,size(conf,2));ones(1,size(conf,2));zeros(1,size(conf,2));zeros(1,size(conf,2));zeros(1,size(conf,2))];
    sensor_circle_config{i} = conf;
end

% Compute the criteria
[reciprocal_condition_number_all_sensor_config, mean_min_singular_value_J_all_sensor_config] = ...
    criteria_from_sensor_config(sensor_circle_config, magnet_conf, B_r, Volume, type);

reciprocal_condition_number_all_sensor_config_3 = reciprocal_condition_number_all_sensor_config;
mean_min_singular_value_J_all_sensor_config_3 = mean_min_singular_value_J_all_sensor_config;



%% Figure 4.10, 4.11
number_of_sensor = [[2:1:10].^2];
number_of_sensor = [number_of_sensor(1),5,number_of_sensor(2:end)];
number_of_sensor_peri = 4:4:100;
number_of_sensor_circle = 4:1:100;

figure
plot(number_of_sensor, reciprocal_condition_number_all_sensor_config_1, '-', 'Color', [0, 0.4470, 0.7410], "linewidth", 2)
hold on
plot(number_of_sensor_peri, reciprocal_condition_number_all_sensor_config_2, '-', 'Color', [0.8500, 0.3250, 0.0980], "linewidth", 2)
hold on 
plot(number_of_sensor_circle, reciprocal_condition_number_all_sensor_config_3, '-', 'Color', [0.9290, 0.6940, 0.1250], "linewidth", 2)
hold on 

xlabel('Number of sensors' , 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
ylabel('$\min_{x\in \mathcal{X}} \frac{1}{\kappa}(x)$', 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
ylim([0,0.3])
grid on

legend('grid 0.32m', 'circle 0.32m', 'circle 0.32m')

figure
plot(number_of_sensor, mean_min_singular_value_J_all_sensor_config_1, '-', 'Color', [0, 0.4470, 0.7410], "linewidth", 2)
hold on
plot(number_of_sensor_peri, mean_min_singular_value_J_all_sensor_config_2, '-', 'Color', [0.8500, 0.3250, 0.0980], "linewidth", 2)
hold on
plot(number_of_sensor_circle, mean_min_singular_value_J_all_sensor_config_3, '-', 'Color', [0.9290, 0.6940, 0.1250], "linewidth", 2)
hold on
legend('grid 0.32m', 'square 0.32m', 'circle 0.32m')

xlabel('Number of sensors', 'FontSize', 12 , 'FontWeight', 'bold' ,'Interpreter', 'latex')
ylabel('$mean_{x \in \mathcal{X}} \sigma_{min}(x)$', 'FontSize', 12 , 'FontWeight', 'bold', 'Interpreter', 'latex')
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
function [reciprocal_condition_number_all_sensor_config, mean_min_singular_value_J_all_sensor_config] = ...
            criteria_from_sensor_config(sensor_config_collection, magnet_conf, B_r, Volume, type)
    % Loop over each sensor configuration
    reciprocal_condition_number_all_sensor_config = [];
    mean_min_singular_value_J_all_sensor_config = [];
    
    for i=1:length(sensor_config_collection)
        
        
        % Extract one sensor configuration
        sensor_config = sensor_config_collection{i};
        sensor_number = size(sensor_config,2);
        sensor_position = sensor_config(1:3,:);
        sensor_orientation = sensor_config(4:7,:);
        magnitudes = vecnorm(sensor_orientation);
        sens_or_unitary = sensor_orientation ./ magnitudes;
    
        % Criteria
        reciprocal_condition_number = [];
        min_singular_value_J = [];
        max_singular_value_J = [];
        mean_singular_value_J = [];
    
        for j=1:size(magnet_conf,2)
            magnet_pos = magnet_conf(1:3,j);
            theta = magnet_conf(4, j);
            phi = magnet_conf(5, j);
            psi = magnet_conf(6, j);
    
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
    
            J = [];
            for k = 1:sensor_number
                J = [J;J_analytical_sensor_B_to_world_pR(sensor_position(:,k), sensor_orientation(:,k), ...
                    magnet_pos, B_r, Volume, R_star, type)];
            end
    
            num_dof = 5;
    
            sigma_J = svd(J);
    
            rcond = sigma_J(num_dof)/sigma_J(1);
            min_svd_J = sigma_J(num_dof);
            max_svd_J = sigma_J(1);
    
            reciprocal_condition_number = [reciprocal_condition_number; rcond];
            min_singular_value_J = [min_singular_value_J; min_svd_J];
    %         max_singular_value_J = [max_singular_value_J; max_svd_J];
            
            % Exclude the 0 sigma
            sigma_J(end) = [];
            mean_singular_value_J = [mean_singular_value_J; mean(sigma_J)];
    
        end
    
        % Put the minimum among all magnet configs into the sensor config
        % collection
        reciprocal_condition_number_all_sensor_config = [reciprocal_condition_number_all_sensor_config; min(reciprocal_condition_number)];
        mean_min_singular_value_J_all_sensor_config = [mean_min_singular_value_J_all_sensor_config; mean(min_singular_value_J)];
    end
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


