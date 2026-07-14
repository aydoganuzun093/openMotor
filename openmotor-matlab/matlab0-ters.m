% Solid Rocket Motor Internal Ballistics Simulator
% BATES (Bore And Tube Ends) Configuration

clear; clc; close all;

%% 1. Input Parameters (Target Goals for Inverse Design)

% Target Performance
target.total_impulse_Ns = 2000;         % Hedeflenen Toplam İtki [Ns]
target.burn_time_s = 2.0;               % Hedeflenen Yanma Süresi [s]
target.chamber_pressure_MPa = 3.5;      % Hedeflenen Ortalama Yanma Odası Basıncı [MPa]

% Input Geometry Constraint
target.outer_diameter_mm = 98.7;        % Sabit Dış Çap [mm]

% Propellant Details (Fixed properties)
propellant.a_burn_rate_coef = 0.025;  % Burn rate coefficient (a) [mm/(s·Pa^n)] (Pa için 0.025)
propellant.n_burn_rate_exp = 0.43;    % Burn rate exponent (n) [-]
propellant.gamma = 1.17;              % Specific heat ratio (k veya gamma) [-]
propellant.Tc_kelvin = 3019;          % Combustion temperature [K]
propellant.mol_mass_kg_kmol = 25.935; % Exhaust molar mass [kg/kmol]
propellant.density_kg_m3 = 1710;      % Propellant density [kg/m^3]

% Nozzle Properties
nozzle.efficiency = 0.95;         % Nozül verimi (Cf düzeltmesi)

% Other Options
other.ambient_pressure_MPa = 0.101325; % Ortam basıncı [MPa]
other.dt_s = 0.001;                    % Zaman adımı [s]

% Assumptions for BATES grain
grain.outer_inhibited = true;     % Dış yüzeyin yanmaz olduğu varsayıldı
grain.ends_inhibited = false;     % Uç kısımlar yanar (BATES şartı)

%% 2. Inverse Design Calculations

% Constants
mm2m = 1e-3;
MPa2Pa = 1e6;
R_u = 8314.46; % Universal gas constant [J/(kmol*K)]
R_spec = R_u / propellant.mol_mass_kg_kmol; % Specific gas constant [J/(kg*K)]
gamma = propellant.gamma;
Tc = propellant.Tc_kelvin;
Pa_atm = other.ambient_pressure_MPa * MPa2Pa;
Pc_target = target.chamber_pressure_MPa * MPa2Pa;

% 2.1 Thermodynamic Properties
% Characteristic Velocity (c*)
c_star = sqrt((gamma * R_spec * Tc) / (gamma^2 * (2/(gamma+1))^((gamma+1)/(gamma-1))));

% Optimum Expansion for Nozzle (Pe = Pa_atm)
% Ideal thrust coefficient (Cf) at optimum expansion
Cf_ideal = sqrt((2*gamma^2/(gamma-1)) * (2/(gamma+1))^((gamma+1)/(gamma-1)) * (1 - (Pa_atm/Pc_target)^((gamma-1)/gamma)));
Cf_delivered = Cf_ideal * nozzle.efficiency;

% Exit Mach Number for Optimum Expansion
Me = sqrt((2/(gamma-1)) * ((Pc_target/Pa_atm)^((gamma-1)/gamma) - 1));

% Optimum Expansion Ratio (Ae / At)
expansion_ratio = (1/Me) * ((2/(gamma+1)) * (1 + (gamma-1)/2 * Me^2))^((gamma+1)/(2*(gamma-1)));

% 2.2 Grain Geometry Calculations
% Burn rate at target chamber pressure
r_burn_mm_s = propellant.a_burn_rate_coef * (Pc_target^propellant.n_burn_rate_exp);

% Required web thickness for the target burn time
web_thickness_mm = r_burn_mm_s * target.burn_time_s;

% Calculate Core Diameter
grain.outer_diameter_mm = target.outer_diameter_mm;
grain.core_diameter_mm = grain.outer_diameter_mm - 2 * web_thickness_mm;

if grain.core_diameter_mm <= 0
    error('Target burn time and pressure result in a core diameter <= 0. Decrease burn time or pressure, or increase outer diameter.');
end

% Single grain length for perfectly neutral burning (BATES Condition: L = 3*r_outer + r_core)
r_outer_mm = grain.outer_diameter_mm / 2;
r_core_mm = grain.core_diameter_mm / 2;
grain.length_mm = 3 * r_outer_mm + r_core_mm;

% Calculate Single Grain Mass
vol_grain_single_mm3 = pi * (r_outer_mm^2 - r_core_mm^2) * grain.length_mm;
vol_grain_single_m3 = vol_grain_single_mm3 * (mm2m^3);
mass_grain_single_kg = vol_grain_single_m3 * propellant.density_kg_m3;

% 2.3 Required Propellant Mass & Number of Grains
% Target Total Mass based on Target Total Impulse
% Itotal = m_prop * c_star * Cf_delivered
target_mass_prop_kg = target.total_impulse_Ns / (c_star * Cf_delivered);

% Calculate Number of Grains (Option 2: Exact Neutrality, Round N)
grain.number = round(target_mass_prop_kg / mass_grain_single_kg);
if grain.number < 1
    grain.number = 1; % At least 1 grain
end

% Actual Propellant Mass with rounded number of grains
actual_mass_prop_kg = grain.number * mass_grain_single_kg;

% 2.4 Nozzle Sizing
% Average Mass Flow Rate based on actual mass and target burn time
mdot_avg = actual_mass_prop_kg / target.burn_time_s;

% Required Throat Area (At = mdot * c* / Pc)
At = (mdot_avg * c_star) / Pc_target;
nozzle.throat_diameter_mm = 2 * sqrt(At / pi) / mm2m;

% Required Exit Area
Ae = At * expansion_ratio;
nozzle.exit_diameter_mm = 2 * sqrt(Ae / pi) / mm2m;

% Print Calculated Design Parameters before Simulation
fprintf('======================================================\n');
fprintf('             INVERSE DESIGN PARAMETERS                \n');
fprintf('======================================================\n');
fprintf('Grain Outer Diameter:      %.2f mm\n', grain.outer_diameter_mm);
fprintf('Grain Core Diameter:       %.2f mm\n', grain.core_diameter_mm);
fprintf('Single Grain Length:       %.2f mm\n', grain.length_mm);
fprintf('Web Thickness:             %.2f mm\n', web_thickness_mm);
fprintf('Number of Grains:          %d\n', grain.number);
fprintf('Nozzle Throat Diameter:    %.2f mm\n', nozzle.throat_diameter_mm);
fprintf('Nozzle Exit Diameter:      %.2f mm\n', nozzle.exit_diameter_mm);
fprintf('Expansion Ratio (Ae/At):   %.2f\n', expansion_ratio);
fprintf('======================================================\n\n');

%% 3. Pre-computations for Simulation
% Grain initial geometry in meters
r_outer = grain.outer_diameter_mm * mm2m / 2;
r_core_init = grain.core_diameter_mm * mm2m / 2;
L_g_init = grain.length_mm * mm2m;
N = grain.number;

% Pressure ratio at exit
P_ratio_exit = (1 + (gamma-1)/2 * Me^2)^(gamma/(gamma-1));

%% 4. Initial Mass and Volume
vol_grain_init = pi * (r_outer^2 - r_core_init^2) * L_g_init;
mass_prop_total = N * vol_grain_init * propellant.density_kg_m3;
chamber_vol_init = N * pi * r_outer^2 * L_g_init; % Simplified total chamber volume
free_vol_init = chamber_vol_init - N * vol_grain_init;

%% 5. Simulation Setup (Iterative Transient Solution)
t_max = 20; % sec, safety limit
num_steps = t_max / other.dt_s;

time_arr = zeros(1, num_steps);
P_arr = zeros(1, num_steps);
F_arr = zeros(1, num_steps);
Kn_arr = zeros(1, num_steps);
mass_flux_arr = zeros(1, num_steps);

% State variables
P = Pa_atm; % Initial pressure is atmospheric (Pa)
free_vol = free_vol_init;
r_core = r_core_init;
L_g = L_g_init;
burned_mass = 0;

P_arr(1) = P;
t = 0;
idx = 1;

ignited = true;

% To avoid infinite loop, burn loop
while ignited && r_core < r_outer && L_g > 0 && idx < num_steps
    % Burn rate (mm/s) expects Pressure in Pa per KNSU (a=0.025 at Pa)
    if P < Pa_atm
        P_burn = Pa_atm;
    else
        P_burn = P;
    end
    r_burn_m_s = (propellant.a_burn_rate_coef * (P_burn^propellant.n_burn_rate_exp)) * mm2m;

    % Current Geometry
    Ab_core = N * 2 * pi * r_core * L_g;
    if grain.ends_inhibited
        Ab_ends = 0;
    else
        Ab_ends = N * 2 * pi * (r_outer^2 - r_core^2); % 2 ends per grain
    end

    Ab_total = Ab_core + Ab_ends;

    % Current Port Area
    Ap = pi * r_core^2;

    % Mass Generation Rate (in)
    mdot_in = propellant.density_kg_m3 * Ab_total * r_burn_m_s;

    % Mass Flow Rate (out)
    if P > Pa_atm
        mdot_out = (P * At) / c_star;
    else
        mdot_out = 0;
    end

    % Transient Pressure Differential Equation
    % dP/dt = (R_spec * Tc / V_c) * (mdot_in - mdot_out) - P * (dV_c / dt) / V_c
    dVc_dt = Ab_total * r_burn_m_s;
    dP_dt = (R_spec * Tc / free_vol) * (mdot_in - mdot_out) - P * dVc_dt / free_vol;

    % Update Pressure (Euler integration)
    P = P + dP_dt * other.dt_s;
    if P < Pa_atm
        P = Pa_atm;
    end

    % Thrust Calculation
    if P > Pa_atm
        Pe = P / P_ratio_exit;
        % Ideal thrust coefficient
        Cf_ideal = sqrt((2*gamma^2/(gamma-1)) * (2/(gamma+1))^((gamma+1)/(gamma-1)) * (1 - (Pe/P)^((gamma-1)/gamma))) + (Pe - Pa_atm)/P * expansion_ratio;
        Cf_delivered = Cf_ideal * nozzle.efficiency;
        Thrust = Cf_delivered * P * At;
    else
        Thrust = 0;
        Cf_ideal = 0;
        Cf_delivered = 0;
    end

    % Update arrays
    time_arr(idx) = t;
    P_arr(idx) = P;
    F_arr(idx) = Thrust;
    Kn_arr(idx) = Ab_total / At;
    mass_flux_arr(idx) = mdot_out / Ap;

    % Update Geometry
    r_core = r_core + r_burn_m_s * other.dt_s;
    if ~grain.ends_inhibited
        L_g = L_g - 2 * r_burn_m_s * other.dt_s;
    end
    free_vol = free_vol + dVc_dt * other.dt_s;

    burned_mass = burned_mass + mdot_in * other.dt_s;

    % Check burnout
    if r_core >= r_outer || L_g <= 0
        ignited = false;
    end

    t = t + other.dt_s;
    idx = idx + 1;
end

% Tail-off (Gas blowdown after burnout)
while P > Pa_atm * 1.05 && idx < num_steps
    mdot_in = 0;
    mdot_out = (P * At) / c_star;
    dVc_dt = 0;
    dP_dt = (R_spec * Tc / free_vol) * (mdot_in - mdot_out) - P * dVc_dt / free_vol;

    P = P + dP_dt * other.dt_s;
    if P < Pa_atm
        P = Pa_atm;
    end

    if P > Pa_atm
        Pe = P / P_ratio_exit;
        Cf_ideal = sqrt((2*gamma^2/(gamma-1)) * (2/(gamma+1))^((gamma+1)/(gamma-1)) * (1 - (Pe/P)^((gamma-1)/gamma))) + (Pe - Pa_atm)/P * expansion_ratio;
        Cf_delivered = Cf_ideal * nozzle.efficiency;
        Thrust = Cf_delivered * P * At;
    else
        Thrust = 0;
    end

    time_arr(idx) = t;
    P_arr(idx) = P;
    F_arr(idx) = Thrust;
    Kn_arr(idx) = 0;
    mass_flux_arr(idx) = mdot_out / Ap;

    t = t + other.dt_s;
    idx = idx + 1;
end

% Trim arrays
time_arr = time_arr(1:idx-1);
P_arr = P_arr(1:idx-1);
F_arr = F_arr(1:idx-1);
Kn_arr = Kn_arr(1:idx-1);
mass_flux_arr = mass_flux_arr(1:idx-1);

%% 6. Calculations & Outputs
% Motor Designation (Total Impulse Classification)
total_impulse = trapz(time_arr, F_arr);
% Motor Designation (Total Impulse Classification)
% Define classes with their lower and upper bounds
motor_classes = {'Micro', 0.0, 0.312; '1/4A', 0.312, 0.625; '1/2A', 0.625, 1.25; ...
                 'A', 1.25, 2.5; 'B', 2.5, 5.0; 'C', 5.0, 10.0; 'D', 10.0, 20.0; ...
                 'E', 20.0, 40.0; 'F', 40.0, 80.0; 'G', 80.0, 160.0; 'H', 160.0, 320.0; ...
                 'I', 320.0, 640.0; 'J', 640.0, 1280.0; 'K', 1280.0, 2560.0; ...
                 'L', 2560.0, 5120.0; 'M', 5120.0, 10240.0; 'N', 10240.0, 20480.0; ...
                 'O', 20480.0, 40960.0; 'P', 40960.0, 81920.0};
motor_letter = 'Unknown';
class_pct = 0;
for i = 1:size(motor_classes, 1)
    low = motor_classes{i, 2};
    high = motor_classes{i, 3};
    if total_impulse > low && total_impulse <= high
        motor_letter = motor_classes{i, 1};
        class_pct = ((total_impulse - low) / (high - low)) * 100;
        break;
    end
end

delivered_isp = total_impulse / (mass_prop_total * 9.80665);

% Burn time: time where Thrust > 5% of Max Thrust
max_thrust = max(F_arr);
thrust_threshold = max_thrust * 0.05;
active_indices = find(F_arr > thrust_threshold);
if ~isempty(active_indices)
    burn_time = time_arr(active_indices(end)) - time_arr(active_indices(1));
else
    burn_time = time_arr(end);
end

volume_loading = (N * vol_grain_init) / chamber_vol_init * 100;
average_pressure_MPa = mean(P_arr(active_indices)) / MPa2Pa;
peak_pressure_MPa = max(P_arr) / MPa2Pa;

initial_Kn = Kn_arr(1);
peak_Kn = max(Kn_arr);

% Ideal thrust coefficient (average during burn)
Pe_avg = mean(P_arr(active_indices)) / P_ratio_exit;
P_avg = mean(P_arr(active_indices));
avg_Cf_ideal = sqrt((2*gamma^2/(gamma-1)) * (2/(gamma+1))^((gamma+1)/(gamma-1)) * (1 - (Pe_avg/P_avg)^((gamma-1)/gamma))) + (Pe_avg - Pa_atm)/P_avg * expansion_ratio;
avg_Cf_delivered = avg_Cf_ideal * nozzle.efficiency;

propellant_length = grain.length_mm * N;
port_throat_ratio = (pi * r_core_init^2) / At;
peak_mass_flux = max(mass_flux_arr);

% Print Results
fprintf('======================================================\n');
fprintf('                SIMULATION RESULTS                    \n');
fprintf('======================================================\n');
fprintf('Motor Designation:         %.1f%% %s\n', class_pct, motor_letter);
fprintf('Impulse:                   %.2f Ns\n', total_impulse);
fprintf('Delivered ISP:             %.2f s\n', delivered_isp);
fprintf('Burn Time:                 %.3f s\n', burn_time);
fprintf('Volume Loading:            %.2f %%\n', volume_loading);
fprintf('Average Pressure:          %.3f MPa\n', average_pressure_MPa);
fprintf('Peak Pressure:             %.3f MPa\n', peak_pressure_MPa);
fprintf('Initial Kn:                %.2f\n', initial_Kn);
fprintf('Peak Kn:                   %.2f\n', peak_Kn);
fprintf('Ideal Thrust Coefficient:  %.3f\n', avg_Cf_ideal);
fprintf('Delivered Thrust Coeff.:   %.3f\n', avg_Cf_delivered);
fprintf('Propellant Mass:           %.3f kg\n', mass_prop_total);
fprintf('Propellant Length:         %.1f mm\n', propellant_length);
fprintf('Port/Throat Ratio:         %.2f\n', port_throat_ratio);
fprintf('Peak Mass Flux:            %.2f kg/(m^2*s)\n', peak_mass_flux);
fprintf('======================================================\n');

%% 7. Plots
figure('Name', 'Thrust and Pressure vs Time', 'NumberTitle', 'off');

% Thrust Plot
subplot(2, 1, 1);
plot(time_arr, F_arr, 'b-', 'LineWidth', 2);
ylabel('Thrust [N]');
xlabel('Time [s]');
title('Motor Performance (Thrust vs. Time)');
ylim([0, max(F_arr)*1.1]);
grid on;

% Pressure Plot
subplot(2, 1, 2);
plot(time_arr, P_arr / MPa2Pa, 'r-', 'LineWidth', 2);
ylabel('Pressure [MPa]');
xlabel('Time [s]');
title('Motor Performance (Pressure vs. Time)');
ylim([0, max(P_arr / MPa2Pa)*1.1]);
grid on;
