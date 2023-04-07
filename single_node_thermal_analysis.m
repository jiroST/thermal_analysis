clear;
clc;

% Satellite's parameters
r_sat = 0.07214;
A_sat = 4 * pi * r_sat^2;
A_sat_sectional = pi * r_sat^2;
m = 1.3;
c_p = 897;
alfa = 0.92;
epsilon = 0.85;

% Satellite's orbit
h_orbit = 550000;

% Constants
stefan_boltzmann_const = 5.6704 * 10^-8;
standard_gravitational_param = 3.986 * 10^14; 
R_earth = 6371000;

% Heat fluxes (W/m^2)
q_sun = 1322;
q_earth = 228;

% Sun's reflection from the earth (%)
albedo = 0.14;

% Rates of heat transfers (W)
Q_insinde = 0.3 * A_sat_sectional * 0.5 * q_sun;
Q_sun = alfa * q_sun * A_sat_sectional;
Q_albedo = albedo * Q_sun;
Q_earth = alfa * q_earth * A_sat_sectional;


% The angle between satellites's orbital plane around Earth and the geocentric position of the sun
beta = 0;

% Beta angle where shifted to eclipse
beta_crit = asind(R_earth / (R_earth + h_orbit));

% Fraction of time in orbit spent in eclipse (%)
if abs(beta) < beta_crit
    eclipse = (1/180) * acosd((sqrt(h_orbit^2 + 2 * R_earth * h_orbit)) / ((R_earth + h_orbit) * cosd(beta)));
elseif abs(beta) >= beta_crit
    eclipse = 0;
end



% Maximum temperature, in the sun
T_max = ((Q_sun + Q_albedo + Q_earth + Q_insinde) / (epsilon * A_sat * stefan_boltzmann_const))^(1/4);
T_max_celsius = T_max - 273.15


% Minimum temperature, in eclipse
T_min = ((Q_earth + Q_insinde) / (epsilon * A_sat * stefan_boltzmann_const))^(1/4);
T_min_celsius = T_min - 273.15


% Average temperature in orbit
T_avg = ((eclipse * (Q_earth + Q_insinde) + (1 - eclipse) * (Q_sun + Q_albedo + Q_earth + Q_insinde)) / ...
        (epsilon * A_sat * stefan_boltzmann_const))^(1/4);
T_avg_celsius = T_avg - 273.15



% Orbital period (s)
tau = 2 * pi * sqrt(((h_orbit + R_earth)^3) / standard_gravitational_param);


% Set up of times and sample sizes for the plotting of the data
n = 3000;
delta_tau = 10*tau/n;
times = linspace(0, 10*tau, n);


% Eclipse timing borders
ecl_border1 = (tau/2) * (1 - eclipse);
ecl_border2 = (tau/2) * (1 + eclipse);

% List for temperatures, let's start from 0°C
temps = [273.15];

for i = 1:1:n-1
    % Switch function turn off solar radiation, when in eclipse 
    switch_func = 0.5 * (sign((cosd((360/tau)*times(i))+0.4))+1);
    
    % Calculated heat transfer in a certain time
    Q_in = q_earth * A_sat_sectional + alfa * (1 + albedo) * q_sun * A_sat_sectional * switch_func;
    Q_radiation = A_sat * epsilon * stefan_boltzmann_const * (temps(i))^4;
    temps(i+1) = temps(i) + (delta_tau / (m * c_p)) * (Q_insinde - Q_radiation + Q_in);
end



% Change units to °C
temps_c = temps - 273.15;


% Plot temperatures
plot(times, temps_c);
hold on;

% Set axis labels and title
xlabel('Time (s)')
ylabel('Temperature (°C)')
title('Single Node Thermal Analysis')

% Draw max, min, and average temperatures
plot([times(1), times(end)], [T_max_celsius, T_max_celsius], 'r--', 'LineWidth', 1.5);
plot([times(1), times(end)], [T_min_celsius, T_min_celsius], 'b--', 'LineWidth', 1.5);
plot([times(1), times(end)], [T_avg_celsius, T_avg_celsius], 'g--', 'LineWidth', 1.5);

% Add legend and grid
legend('Temperature', 'Max temperature', 'Min temperature', 'Avg temperature');
grid on;







