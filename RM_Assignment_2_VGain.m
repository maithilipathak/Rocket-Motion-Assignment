% Assignment 2 : ΔV-GAIN ANALYSIS FOR AIR-LAUNCH SCENARIOS

clear all; close all; clc;

%% Physical Constants
fprintf('ΔV-GAIN ANALYSIS FOR AIR-LAUNCH SCENARIOS\n');


% Earth parameters
R_earth = 6371*1000;          % Earth radius [m]
mu = 3.986004418e14;      % Earth gravitational parameter [m³/s²]


% Atmosphere parameters (US Standard Atmosphere 1976)
T0 = 288.15;             % Sea level temperature [K]
L = 0.0065;              % Temperature lapse rate [K/m] in troposphere
R_air = 8.314;          % Specific gas constant for air [J/(mol·K)]
gamma_air = 1.4;         % Ratio of specific heats for air
M0 = 0.02897;              %Molar mass of air [kg/mol]

% Target orbit
h_target = 185*1000;          % Target orbit altitude [m]
r_target = R_earth + h_target;  % Target orbital radius [m]

%% Define Initial Conditions
% Initial altitudes [km]
h0_values = [0, 7.5, 15, 22.5];

% Initial Mach numbers
Mach_values = [0, 1, 2, 3];

% Number of cases
n_altitudes = length(h0_values);
n_mach = length(Mach_values);
n_cases = n_altitudes * n_mach;

%% Step 1: Calculate Atmospheric Properties and Speed of Sound

fprintf('ATMOSPHERIC PROPERTIES\n');


% Preallocate arrays
T_h0 = zeros(size(h0_values));      % Temperature at altitude [K]
a_h0 = zeros(size(h0_values));      % Speed of sound at altitude [m/s]

fprintf('Altitude [km] | Temperature [K] | Speed of Sound [m/s] \n');
fprintf('--------------------------------------------------------\n');

for i = 1:n_altitudes
    h = h0_values(i);  % Altitude in km
    
    % Temperature at altitude (valid for troposphere, h < 11 km)
    % T(h) = T₀ - L·h
    if h < 11
        T_h0(i) = T0 - (L * h * 1000);  % Convert km to m for lapse rate
    else
        % Stratosphere (isothermal at 11-25 km)
        T_h0(i) = T0 - (L * 11000);  % Temperature at tropopause
    end
    
    % Speed of sound: a = sqrt(γ·R·T/M)
    a_h0(i) = sqrt((gamma_air * R_air * T_h0(i))/M0);  % [m/s]
    
    fprintf('   %5.1f      |     %6.2f      |       %6.2f         \n', ...
            h, T_h0(i), a_h0(i));
end

%% Step 2: Calculate Initial Velocities

fprintf('INITIAL VELOCITIES (V₀ = M₀ × a)\n');

% Create velocity matrix [n_mach × n_altitudes]
V0_matrix = zeros(n_mach, n_altitudes);

fprintf('         | h₀ = 0 km  | h₀ = 7.5 km | h₀ = 15 km  | h₀ = 22.5 km\n');
fprintf('----------------------------------------------------------------------\n');

for i = 1:n_mach
    M = Mach_values(i);
    fprintf('M₀ = %d   |', M);
    
    for j = 1:n_altitudes
        V0_matrix(i,j) = M * a_h0(j);  % V₀ = M₀ × a [m/s]
        fprintf('  %.4f km/s |', V0_matrix(i,j));
    end
    fprintf('\n');
end
%all values are correct till here. cross checked with https://www.engineeringtoolbox.com/elevation-speed-sound-air-d_1534.html
%% Step 3: Calculate Specific Orbital Energy



%taking Specific energy of ground launch as reference.
epsilon_ground = -mu /R_earth;  % Specific energy [km²/s²]
fprintf('Reference Ground Launch Specific Energy (m^2/s^2) %d',epsilon_ground);

% Initial specific energies for all cases
%epsilon = V0^2 / 2 - mu / r0
epsilon_target = zeros(n_mach, n_altitudes);

fprintf('\n');
fprintf('E0 Specific Energies:\n');
fprintf('         | h₀ = 0 km    | h₀ = 7.5 km   | h₀ = 15 km    | h₀ = 22.5 km\n');
fprintf('----------------------------------------------------------------------\n');

for i = 1:n_mach
    M = Mach_values(i);
    fprintf('M₀ = %d   |', M);
    
    for j = 1:n_altitudes
        h0 = h0_values(j);
        r0 = R_earth + h0*1000;  % Initial radius [m]
        V0 = V0_matrix(i,j);
        
        % Specific energy: kinetic + potential
        epsilon_target(i,j) = V0^2 / 2 - mu / (r0);  % [m²/s²]
        
        fprintf('  %8.4f    |', epsilon_target(i,j));
    end
    fprintf('\n');
end

%% Step 4: Energy Difference Required

fprintf('ENERGY DIFFERENCE REQUIRED (Δε = ε_target - ε_initial)\n');


Delta_epsilon = zeros(n_mach, n_altitudes);

fprintf('         | h₀ = 0 km  | h₀ = 7.5 km | h₀ = 15 km  | h₀ = 22.5 km\n');
fprintf('----------------------------------------------------------------------\n');

for i = 1:n_mach
    M = Mach_values(i);
    fprintf('M₀ = %d   |', M);
    
    for j = 1:n_altitudes
        Delta_epsilon(i,j) = epsilon_target(i,j) - epsilon_ground;  % [m²/s²]
        fprintf('  %8.4f   |', Delta_epsilon(i,j));
    end
    fprintf('\n');
end

%% Step 5: Calculate Equivalent ΔV

fprintf('EQUIVALENT ΔV CALCULATION\n');

Delta_V = zeros(n_mach, n_altitudes);

fprintf('Equivalent ΔV Required:\n');
fprintf('         | h₀ = 0 km  | h₀ = 7.5 km | h₀ = 15 km  | h₀ = 22.5 km\n');
fprintf('----------------------------------------------------------------------\n');

for i = 1:n_mach
    M = Mach_values(i);
    fprintf('M₀ = %d   |', M);
    
    for j = 1:n_altitudes
        v_i = V0_matrix(i,j);  % Initial velocity
        Delta_eps = Delta_epsilon(i,j);
        
        % Calculate required final velocity to achieve target energy
             % v_f = √(2·ε_target )
        Delta_V(i,j) = sqrt(2 * Delta_eps);
        
              
        fprintf('  %7.4f  |', Delta_V(i,j));
    end
    fprintf('\n');
end
