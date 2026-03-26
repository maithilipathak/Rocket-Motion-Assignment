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

%% Calculate Atmospheric Properties and Speed of Sound

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

%% Calculate Initial Velocities

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


%% TOTAL SPECIFIC ENERGY AND ΔV CALCULATION


% Target orbit energy (circular orbit)
v_target_orbit = sqrt(mu / r_target);  % Circular velocity at target
epsilon_target_orbit = v_target_orbit^2 / 2 - mu / r_target;  % Total specific energy

% Initial specific energies for each altitude (starting from rest, v₀=0)
epsilon_initial = zeros(1, n_altitudes);



for j = 1:n_altitudes
    h0 = h0_values(j);
    r0 = R_earth + h0*1000;
    % ε = v²/2 - μ/r, with v=0 for altitude-only contribution
    epsilon_initial(j) = 0^2 / 2 - mu / r0;
end

% Energy differences needed to reach target orbit
Delta_epsilon = epsilon_target_orbit - epsilon_initial;  % [1×4] vector

% Convert energy differences to required ΔV
% Using: ΔV = √(v₀² + 2Δε) - v₀
% For v₀=0: ΔV = √(2Δε)
DV_required = sqrt(2 * Delta_epsilon);  % [1×4] vector [m/s]

%% Step 5: Calculate Equivalent ΔV (Potential + Kinetic Energy)

% Reference: ΔV required from ground launch
DV_reference = DV_required(1);  % Ground launch (h=0)

% ΔV-gain from altitude = DV_ground - DV_altitude
DV_gain_altitude = DV_reference - DV_required;  % [1×4] vector [m/s]

%% Step 6: ADD VELOCITY CONTRIBUTION SEPARATELY

% The velocity contribution is simply the initial velocity itself
% DV_gain_velocity = V₀
DV_gain_velocity = V0_matrix;  % [4×4] matrix [m/s]

%% Step 7: TOTAL ΔV-GAIN (Sum of Both Contributions)

% Total ΔV-gain: altitude contribution + velocity contribution
% DV_gain_altitude is [1×4], will broadcast across all rows
% DV_gain_velocity is [4×4]
Delta_V_gain_total = DV_gain_altitude + DV_gain_velocity;  % [4×4] matrix [m/s]

%% Step 8: FINAL RESULTS TABLE

fprintf('TOTAL ΔV-GAIN\n');


fprintf('Altitude [km] | Mach | Velocity [m/s] | Delta V gain [m/s]\n');
fprintf('---------------------------------------------------------------\n');

for j = 1:n_altitudes
    for i = 1:n_mach
        fprintf('    %5.1f     |  %d   |    %10.2f  |    %14.6f\n', ...
                h0_values(j), Mach_values(i), V0_matrix(i,j), ...
                Delta_V_gain_total(i,j));
    end
end
