% =========================================================================
% Assignment 2 : ΔV-GAIN ANALYSIS FOR AIR-LAUNCH SCENARIOS

clear all; close all; clc;

%% Physical Constants
fprintf('ΔV-GAIN ANALYSIS FOR AIR-LAUNCH SCENARIOS\n');


% Earth parameters
R_earth = 6371;          % Earth radius [km]
mu = 3.986004418e5;      % Earth gravitational parameter [km³/s²]
g0 = 9.80665e-3;         % Sea-level gravity [km/s²]

% Atmosphere parameters (US Standard Atmosphere 1976)
T0 = 288.15;             % Sea level temperature [K]
L = 0.0065;              % Temperature lapse rate [K/m] in troposphere
R_gas = 287.05;          % Specific gas constant for air [J/(kg·K)]
gamma_air = 1.4;         % Ratio of specific heats for air

% Target orbit
h_target = 185;          % Target orbit altitude [km]
r_target = R_earth + h_target;  % Target orbital radius [km]

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
a_h0 = zeros(size(h0_values));      % Speed of sound at altitude [km/s]

fprintf('Altitude [km] | Temperature [K] | Speed of Sound [m/s] \n');
fprintf('--------------------------------------------------------\n');

for i = 1:n_altitudes
    h = h0_values(i);  % Altitude in km
    
    % Temperature at altitude (valid for troposphere, h < 11 km)
    % T(h) = T₀ - L·h
    if h <= 11
        T_h0(i) = T0 - L * h * 1000;  % Convert km to m for lapse rate
    else
        % Stratosphere (isothermal at 11-25 km)
        T_h0(i) = T0 - L * 11000;  % Temperature at tropopause
    end
    
    % Speed of sound: a = sqrt(γ·R·T)
    a_h0(i) = sqrt(gamma_air * R_gas * T_h0(i)) / 1000;  % [km/s]
    
    fprintf('   %5.1f      |     %6.2f      |       %6.2f         \n', ...
            h, T_h0(i), a_h0(i)*1000);
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
        V0_matrix(i,j) = M * a_h0(j);  % V₀ = M₀ × a [km/s]
        fprintf('  %.4f km/s |', V0_matrix(i,j));
    end
    fprintf('\n');
end

%% Step 3: Calculate Specific Orbital Energy

fprintf('SPECIFIC ORBITAL ENERGY\n');

% Target orbit specific energy (circular orbit)
% For circular orbit: v_circular = sqrt(μ/r)
v_target = sqrt(mu / r_target);  % Circular orbital velocity [km/s]
epsilon_target = v_target^2 / 2 - mu / r_target;  % Specific energy [km²/s²]

fprintf('Target Orbit (circular at h = %d km):\n', h_target);
fprintf('  Orbital radius: r = %.0f km\n', r_target);
fprintf('  Circular velocity: v_circ = √(μ/r) = %.4f km/s = %.0f m/s\n', ...
        v_target, v_target*1000);
fprintf('  Specific energy: ε_target = v²/2 - μ/r = %.4f km²/s²\n\n', epsilon_target);

% Initial specific energies for all cases
epsilon_initial = zeros(n_mach, n_altitudes);
r_initial = zeros(1, n_altitudes);

fprintf('Initial Specific Energies:\n');
fprintf('         | h₀ = 0 km  | h₀ = 7.5 km | h₀ = 15 km  | h₀ = 22.5 km\n');
fprintf('----------------------------------------------------------------------\n');

for i = 1:n_mach
    M = Mach_values(i);
    fprintf('M₀ = %d   |', M);
    
    for j = 1:n_altitudes
        h0 = h0_values(j);
        r0 = R_earth + h0;  % Initial radius [km]
        r_initial(j) = r0;
        
        V0 = V0_matrix(i,j);
        
        % Specific energy: kinetic + potential
        epsilon_initial(i,j) = V0^2 / 2 - mu / r0;  % [km²/s²]
        
        fprintf('  %8.4f   |', epsilon_initial(i,j));
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
        Delta_epsilon(i,j) = epsilon_target - epsilon_initial(i,j);  % [km²/s²]
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
        % ε_target = v_f²/2 - μ/r_i
        % v_f = √(2·ε_target + 2μ/r_i) = √(2·Δε + v_i²)
        v_f = sqrt(2 * Delta_eps + v_i^2);
        
        % ΔV is the difference
        Delta_V(i,j) = v_f - v_i;  % [km/s]
        
        fprintf('  %7.4f  |', Delta_V(i,j));
    end
    fprintf('\n');
end

%% Step 6: Calculate ΔV-Gain Relative to Ground Launch

fprintf('ΔV-GAIN RELATIVE TO GROUND LAUNCH AT REST\n');


% Reference case: ground launch at rest (M=0, h=0)
Delta_V_reference = Delta_V(1,1);  % M₀=0, h₀=0 km

fprintf('Reference Case (Ground Launch at Rest):\n');
fprintf('  h₀ = 0 km, M₀ = 0\n');
fprintf('  Required ΔV = %.4f km/s = %.0f m/s\n\n', Delta_V_reference, Delta_V_reference*1000);

% Calculate ΔV-gain
Delta_V_gain = zeros(n_mach, n_altitudes);

fprintf('ΔV-GAIN = ΔV_reference - ΔV_case [km/s]:\n');
fprintf('         | h₀ = 0 km  | h₀ = 7.5 km | h₀ = 15 km  | h₀ = 22.5 km\n');
fprintf('----------------------------------------------------------------------\n');

for i = 1:n_mach
    M = Mach_values(i);
    fprintf('M₀ = %d   |', M);
    
    for j = 1:n_altitudes
        Delta_V_gain(i,j) = Delta_V_reference - Delta_V(i,j);  % [km/s]
        fprintf('  %7.4f  |', Delta_V_gain(i,j));
    end
    fprintf('\n');
end

fprintf('\nΔV-GAIN as Percentage of Reference:\n');
fprintf('         | h₀ = 0 km  | h₀ = 7.5 km | h₀ = 15 km  | h₀ = 22.5 km\n');
fprintf('----------------------------------------------------------------------\n');

for i = 1:n_mach
    M = Mach_values(i);
    fprintf('M₀ = %d   |', M);
    
    for j = 1:n_altitudes
        gain_percent = 100 * Delta_V_gain(i,j) / Delta_V_reference;
        fprintf('   %6.2f%% |', gain_percent);
    end
    fprintf('\n');
end

fprintf('\n');

%% Step 7: Detailed Analysis and Visualization

fprintf('DETAILED ANALYSIS\n');

% Find best and worst cases
[max_gain, max_idx] = max(Delta_V_gain(:));
[max_i, max_j] = ind2sub(size(Delta_V_gain), max_idx);

fprintf('MAXIMUM ΔV-GAIN:\n');
fprintf('  Case: M₀ = %d, h₀ = %.1f km\n', Mach_values(max_i), h0_values(max_j));
fprintf('  Initial velocity: V₀ = %.4f km/s = %.0f m/s\n', ...
        V0_matrix(max_i, max_j), V0_matrix(max_i, max_j)*1000);
fprintf('  ΔV required: %.4f km/s\n', Delta_V(max_i, max_j));
fprintf('  ΔV-gain: %.4f km/s = %.0f m/s\n', max_gain, max_gain*1000);
fprintf('  Percentage gain: %.2f%%\n\n', 100*max_gain/Delta_V_reference);

% Analyze contributions
fprintf('BREAKDOWN OF ΔV-GAIN SOURCES:\n\n');

% Pure altitude effect (M=0)
fprintf('1. Altitude Effect Only (M₀ = 0):\n');
for j = 2:n_altitudes
    h = h0_values(j);
    gain = Delta_V_gain(1,j);
    fprintf('   h₀ = %5.1f km: ΔV-gain = %.4f km/s (%.2f%%)\n', ...
            h, gain, 100*gain/Delta_V_reference);
end
fprintf('   Conclusion: Altitude alone provides %.2f%% gain at 22.5 km\n\n', ...
        100*Delta_V_gain(1,4)/Delta_V_reference);

% Pure velocity effect (h=0)
fprintf('2. Velocity Effect Only (h₀ = 0 km):\n');
for i = 2:n_mach
    M = Mach_values(i);
    gain = Delta_V_gain(i,1);
    fprintf('   M₀ = %d: ΔV-gain = %.4f km/s (%.2f%%)\n', ...
            M, gain, 100*gain/Delta_V_reference);
end
fprintf('   Conclusion: Mach 3 at sea level provides %.2f%% gain\n\n', ...
        100*Delta_V_gain(4,1)/Delta_V_reference);

% Combined effect
fprintf('3. Combined Effect (M₀ = 3, h₀ = 22.5 km):\n');
gain_combined = Delta_V_gain(4,4);
gain_alt = Delta_V_gain(1,4);
gain_vel = Delta_V_gain(4,1);
fprintf('   Total gain: %.4f km/s (%.2f%%)\n', gain_combined, 100*gain_combined/Delta_V_reference);
fprintf('   Altitude contribution: %.4f km/s\n', gain_alt);
fprintf('   Velocity contribution: %.4f km/s\n', gain_vel);
fprintf('   Synergy: %.4f km/s\n', gain_combined - gain_alt - gain_vel);


%% Visualization

% Figure 1: Individual plots for each Mach number (2x2 grid)
figure('Position', [50, 50, 1200, 900], 'Color', 'w');

colors_altitude = [0 0.4470 0.7410;    % Blue
                   0.8500 0.3250 0.0980;  % Red-orange
                   0.9290 0.6940 0.1250;  % Yellow-orange
                   0.4940 0.1840 0.5560]; % Purple

for i = 1:n_mach
    subplot(2, 2, i);
    hold on; grid on; box on;
    
    M = Mach_values(i);
    
    % Plot ΔV-gain vs altitude for this Mach number
    plot(h0_values, Delta_V_gain(i,:)*1000, 'o-', ...
         'LineWidth', 2.5, 'MarkerSize', 10, ...
         'Color', colors_altitude(i,:), ...
         'MarkerFaceColor', colors_altitude(i,:), ...
         'MarkerEdgeColor', 'k');
    
    % Add value labels on each point
    for j = 1:n_altitudes
        text(h0_values(j), Delta_V_gain(i,j)*1000 + 25, ...
             sprintf('%.0f m/s', Delta_V_gain(i,j)*1000), ...
             'HorizontalAlignment', 'center', 'FontSize', 10, ...
             'FontWeight', 'bold', 'Color', colors_altitude(i,:));
    end
    
    % Formatting
    xlabel('Initial Altitude h_0 [km]', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('\DeltaV-Gain [m/s]', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Mach %d (V_0 = %.0f - %.0f m/s)', M, ...
          V0_matrix(i,1)*1000, V0_matrix(i,end)*1000), ...
          'FontSize', 13, 'FontWeight', 'bold');
    
    % Set consistent y-axis limits for comparison
    ylim([0, 1100]);
    xlim([-2, 24]);
    set(gca, 'FontSize', 11, 'LineWidth', 1.2);
    
    % Add reference line at y=0
    plot([h0_values(1), h0_values(end)], [0, 0], 'k--', 'LineWidth', 1);
    
    % Add percentage gain annotation
    gain_percent = 100 * Delta_V_gain(i,:) / Delta_V_reference;
    text(h0_values(end)-2, Delta_V_gain(i,end)*1000 - 80, ...
         sprintf('Max: %.1f%%', gain_percent(end)), ...
         'FontSize', 11, 'FontWeight', 'bold', ...
         'Color', [0.5 0.5 0.5], ...
         'HorizontalAlignment', 'right');
end

sgtitle('\DeltaV-Gain vs Altitude for Different Mach Numbers', ...
        'FontSize', 16, 'FontWeight', 'bold');

%% Figure 2: Combined comparison plot
figure('Position', [100, 100, 1000, 700], 'Color', 'w');

subplot(2,1,1);
hold on; grid on; box on;

% Plot all Mach numbers on same graph
markers = {'o', 's', '^', 'd'};
colors = [0 0.4470 0.7410;    % Blue (Mach 0)
          0.8500 0.3250 0.0980;  % Red (Mach 1)
          0.9290 0.6940 0.1250;  % Orange (Mach 2)
          0.4940 0.1840 0.5560]; % Purple (Mach 3)

for i = 1:n_mach
    plot(h0_values, Delta_V_gain(i,:)*1000, [markers{i} '-'], ...
         'LineWidth', 2.5, 'MarkerSize', 10, ...
         'Color', colors(i,:), ...
         'MarkerFaceColor', colors(i,:), ...
         'MarkerEdgeColor', 'k', ...
         'DisplayName', sprintf('Mach %d', Mach_values(i)));
end

xlabel('Initial Altitude h_0 [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('\DeltaV-Gain [m/s]', 'FontSize', 13, 'FontWeight', 'bold');
title('\DeltaV-Gain vs Altitude (All Mach Numbers)', ...
      'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 11);
ylim([0, 1100]);
set(gca, 'FontSize', 12, 'LineWidth', 1.2);

% Subplot 2: ΔV-Gain as percentage
subplot(2,1,2);
hold on; grid on; box on;

gain_percent_matrix = 100 * Delta_V_gain / Delta_V_reference;

for i = 1:n_mach
    plot(h0_values, gain_percent_matrix(i,:), [markers{i} '-'], ...
         'LineWidth', 2.5, 'MarkerSize', 10, ...
         'Color', colors(i,:), ...
         'MarkerFaceColor', colors(i,:), ...
         'MarkerEdgeColor', 'k', ...
         'DisplayName', sprintf('Mach %d', Mach_values(i)));
end

xlabel('Initial Altitude h_0 [km]', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('\DeltaV-Gain [% of Reference]', 'FontSize', 13, 'FontWeight', 'bold');
title('\DeltaV-Gain as Percentage of Ground Launch', ...
      'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 11);
ylim([0, 14]);
set(gca, 'FontSize', 12, 'LineWidth', 1.2);

%% Summary Table

fprintf('FINAL SUMMARY TABLE\n');

fprintf('┌────────┬─────────┬────────────┬──────────────┬──────────────┬────────────┐\n');
fprintf('│  Case  │   h₀    │     M₀     │   V₀ [m/s]   │  ΔV [m/s]    │ Gain [m/s] │\n');
fprintf('├────────┼─────────┼────────────┼──────────────┼──────────────┼────────────┤\n');

case_num = 1;
for j = 1:n_altitudes
    for i = 1:n_mach
        fprintf('│  %2d    │ %5.1f   │    %d       │   %7.1f    │   %7.1f    │  %7.1f   │\n', ...
                case_num, h0_values(j), Mach_values(i), ...
                V0_matrix(i,j)*1000, Delta_V(i,j)*1000, Delta_V_gain(i,j)*1000);
        case_num = case_num + 1;
    end
    if j < n_altitudes
        fprintf('├────────┼─────────┼────────────┼──────────────┼──────────────┼────────────┤\n');
    end
end
fprintf('└────────┴─────────┴────────────┴──────────────┴──────────────┴────────────┘\n\n');

fprintf('Reference (Ground Launch at Rest): ΔV = %.0f m/s\n', Delta_V_reference*1000);
fprintf('Maximum Gain: %.0f m/s (%.1f%%) at M₀=%d, h₀=%.1f km\n\n', ...
        max_gain*1000, 100*max_gain/Delta_V_reference, ...
        Mach_values(max_i), h0_values(max_j));
