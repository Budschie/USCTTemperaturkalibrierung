% Obtained from https://iopscience.iop.org/article/10.1088/1742-6596/279/1/012023/pdf
graphics_toolkit("fltk");

% expected_weight_perc = 12.6;
expected_weight_perc = 15.595;
expected_weight_perc = 15;
% set_temperature = 36.2;
set_temperature = 20;

global min_temp max_temp;

%min_temp = 35.2 - 4;
% max_temp = 37.4 + 4;

% Avg body temperature with mean age of 30
average_women_body_temperature = 36.2;

min_temp = 35.2;
max_temp = 37.2;

% min_temp = 18;
% max_temp = 45;

% min_temp = average_women_body_temperature - 1;
% max_temp = average_women_body_temperature + 1;
% min_temp = 19;
% max_temp = 21;

% min_temp = 18;
% max_temp = 28;
% temp_resolution = 0.001;
temp_resolution = 0.03;

calculate_functions

function c_list = get_c_for_temperature_range_and_wp(temperature_range, wp)
  # Initialize vector that can hold all c values
  c_list = zeros(size(temperature_range)(2), 1);

  # WHY ON EARTH DOES THIS HAVE INDEX 1?
  size_test = size(c_list)(1);

  for i = 1:size_test
    c_list(i) = calculate_speed_of_sound(temperature_range(i), wp);
  endfor
end

function [mu, sigma] = calc_default_statistics(c_list)
  [z, mu, sigma] = zscore(c_list);
endfunction

function [mu, sigma] = calc_default_statistics_args(temperature_range, wp)
  [mu, sigma] = calc_default_statistics(get_c_for_temperature_range_and_wp(temperature_range, wp));
endfunction

changing_temperature = min_temp:temp_resolution:max_temp;
% sos_for_temp = zeros(size(changing_temperature)(2), 1);

% for i = 1:size(sos_for_temp)(1)
%  sos_for_temp(i) = calculate_speed_of_sound(changing_temperature(i), expected_weight_perc);
% endfor

sos_for_temp = calculate_speed_of_sound(changing_temperature, repmat(expected_weight_perc, length(changing_temperature), 1));

% changing_wp = expected_weight_perc-4:0.002:expected_weight_perc+4;
changing_wp = expected_weight_perc-4:0.02:expected_weight_perc+4;
diff_sos = zeros(size(changing_wp)(2), 1);

set_c = calculate_speed_of_sound(set_temperature, expected_weight_perc);

for i = 1:size(diff_sos)(1)
  diff_sos(i) = calculate_speed_of_sound(set_temperature, changing_wp(i)) - set_c;
endfor

% changing_wp_min_max_finder = 0:0.005:25;
changing_wp_min_max_finder = 0:0.0005:25;
temperature_dependence_metric = zeros(size(changing_wp_min_max_finder)(2), 1);
% sigma_values = zeros(size(changing_wp_min_max_finder)(2), 1);
% mu_values = zeros(size(changing_wp_min_max_finder)(2), 1);

function abs_diff = get_mean_abs_diff(diff_values, mu)
  abs_diff = mean(abs(diff_values - mu));
endfunction

function [metric, mu] = get_temperature_dependence_metric(temperature_range, wp)
  % all_c_values = calculate_bulk_speed_of_sound(temperature_range, repmat(wp, length(temperature_range), 1));
  [temperature_range_grid, wp_grid] = meshgrid(temperature_range, wp);
  all_c_values = calculate_speed_of_sound(temperature_range_grid, wp_grid)';

  median_val = get_mean_abs_diff(all_c_values, median(all_c_values));
  mean_val = get_mean_abs_diff(all_c_values, mean(all_c_values));

  trial_error_val = 10000;

  % for i = 1:length(all_c_values)
  % trial_error_val = min(trial_error_val, get_mean_abs_diff(all_c_values, all_c_values(i)));
  % endfor

  smallest = min(min(median_val, mean_val), trial_error_val);

  %if smallest == trial_error_val
  %  printf("Trial method\n");
  %endif

  %if smallest == mean_val
  %  printf("Mean method\n");
  %endif

  % if smallest == median_val
  %  printf("Median method\n");
  %endif

  % if trial_error_val == median_val
  %  printf("YES\n");
  % else
  %  printf("NOPE\n");
  % endif

  % printf("MEAN: %.4e; MEDIAN: %.4e", mean_val, median_val);

  % Median seems to be the most optimal value, so it is used as a metric
  metric = smallest;
  mu = median(all_c_values);
endfunction

[temperature_dependence_metric, mu_tdm] = get_temperature_dependence_metric(changing_temperature, changing_wp_min_max_finder);

% for i = 1:length(temperature_dependence_metric)
  % temperature_dependence_metric(i) = get_temperature_dependence_metric(changing_temperature, changing_wp_min_max_finder(i));
  % all_c_values = get_c_for_temperature_range_and_wp(changing_temperature, changing_wp_min_max_finder(i));
  % max_val = max(all_c_values);
  % min_val = min(all_c_values);
  % [mu, sigma] = calc_default_statistics_args(changing_temperature, changing_wp_min_max_finder(i));
  % sigma_values(i) = sigma;
  % mu_values(i) = mu;
% endfor

wp_range = 0:0.001:25;
c_against_wp = zeros(length(wp_range), 1);

for i = 1:length(wp_range)
  c_against_wp(i) = calculate_speed_of_sound(set_temperature, wp_range(i));
endfor

c_against_wp = calculate_speed_of_sound(set_temperature, wp_range);

[min_y, min_x] = min(temperature_dependence_metric);

baseline_c_values = calculate_bulk_speed_of_sound(changing_temperature, repmat(expected_weight_perc, length(changing_temperature), 1));

function standard_error = get_standard_error(expected, estimate)
  mse_single = (expected - estimate) .* (expected - estimate);
  mse_averaged = mean(mse_single, 1);

  standard_error = sqrt(mse_averaged);
endfunction

% xi_deviations = -4:0.001:4;
xi_deviations = -4:0.1:4;

error_metrics = zeros(length(xi_deviations), 1);

for i = 1:length(error_metrics)
  error_metrics(i) = get_standard_error(baseline_c_values, calculate_bulk_speed_of_sound(changing_temperature, repmat(expected_weight_perc + xi_deviations(i), length(changing_temperature), 1)));
  % error_metrics(i) = get_temperature_dependence_metric(temperature_range, expected_weight_perc + xi_deviations(i));
endfor

wp_used = changing_wp_min_max_finder(min_x);
printf("c in w%%/w%% %.6f; σ in m/s %.8e µ_opt=%.8e\n", wp_used, min_y, mu_tdm(min_x));

printf("pure water; σ in m/s %.8e µ_opt=%.8e\n", temperature_dependence_metric(1), mu_tdm(1));

printf("=== STD DEV STATS ===\n");

[mu_wp, sigma_wp] = calc_default_statistics_args(min_temp:temp_resolution:max_temp, wp_used);
printf("c for w%%/w%% %.3f has statistical properties squared error µ=%.8e and σ=%.8e \n", wp_used, mu_wp, sigma_wp);

[mu_water, sigma_water] = calc_default_statistics_args(min_temp:temp_resolution:max_temp, 0);
printf("c for pure water has statistical properties squared error µ=%.8e and σ=%.8e \n", mu_water, sigma_water);

printf("=== END ===\n");
printf("=== END ===\n");

figure(1);
plot(changing_temperature, sos_for_temp, "linewidth", 8);
set(gca, "linewidth", 6, "fontsize", 22);
title(sprintf("Speed of sound against temperature for 𝛏=%.2f%%", expected_weight_perc));
xlabel("T in °C");
ylabel("c in m/s");

figure(2);
plot(changing_wp, diff_sos, "linewidth", 8);
set(gca, "linewidth", 6, "fontsize", 22);
title(sprintf("Difference between predicted speed of sound at wp=%.3f%% and real wp at T=%.3f°C", expected_weight_perc, set_temperature));
xlabel("WP");
ylabel("c in m/s");

figure(3);
plot(changing_wp_min_max_finder, temperature_dependence_metric, "linewidth", 8);
set(gca, "linewidth", 6, "fontsize", 22);
title(sprintf("σ von c mit %.2f°C ≤ T ≤ %.2f°C gegen 𝛏", min_temp, max_temp));
xlabel("𝛏 in Prozent");
ylabel("σ in m/s");

figure(4);
plot(wp_range, c_against_wp, "linewidth", 8);
set(gca, "linewidth", 6, "fontsize", 22);
title(sprintf("Change of c in regards to changing wp at constant T=%.3f°C", set_temperature));
xlabel("WP");
ylabel("c in m/s");

figure(5);
plot(xi_deviations, error_metrics, "linewidth", 3);

hold on;
plot(xi_deviations, repmat(1.0, length(xi_deviations), 1), "linewidth", 1);
plot(xi_deviations, repmat(0.1, length(xi_deviations), 1), "linewidth", 1);
hold off;
% set(gca, "linewidth", 6, "fontsize", 22);
title(sprintf("Standard error at assumed $\xi$=%.3f%% due to wrong $\xi$", expected_weight_perc));
xlabel("$\Delta \xi$ in absolute weight percentage");
ylabel("Standard error of expected vs real $c$ in m/s");

legend("Standard error", "Minimum required accuracy", "1/10 of required accuracy");
% matlab2tikz("wp_error.tex", "showInfo", false, 'height', '\fheight', 'width', '\fwidth');

