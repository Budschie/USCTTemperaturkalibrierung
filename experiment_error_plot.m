% Obtained from https://iopscience.iop.org/article/10.1088/1742-6596/279/1/012023/pdf
graphics_toolkit("qt");

global DI_INDEX KI_INDEX LI_INDEX speed_of_sound_matrix;

DI_INDEX = 2;
KI_INDEX = 3;
LI_INDEX = 4;

b = 0.015479;
vk = 2552;
% Rounding error
error = 1e-8;

set_temperature = 30;
set_weight_percentage = 13.9;
temp_resolution = 0.05;

speed_of_sound_matrix = dlmread("UltrasonicSpeedList.csv", ",");

function speed_of_sound = calculate_speed_of_sound(temperature, ethanol_weight_percentage)
  global DI_INDEX KI_INDEX LI_INDEX speed_of_sound_matrix;
  e_power_k_i = ethanol_weight_percentage .^ speed_of_sound_matrix(:,KI_INDEX);
  t_power_l_i = temperature .^ speed_of_sound_matrix(:,LI_INDEX);

  catted = cat(2, e_power_k_i, t_power_l_i, speed_of_sound_matrix(:,DI_INDEX));

  speed_of_sound_elements = prod(catted, 2);
  speed_of_sound = sum(speed_of_sound_elements);
end

% Errors resulting from the rounding error above. Naming is bad, ik...
function minus_error = calculate_negative_error(error, value)
  minus_error = -(error) / (value * (value - error));
end

function plus_error = calculate_positive_error(error, value)
  plus_error = (error) / (value * (value + error));
end

function error_skew = calculate_error_skew(error, value)
  error_skew = abs(value - error) / (error + value);
end

function value = get_value_v_known(b, v)
  value = b / v;
end

function value = get_value_time_difference_known(b, vk, time_diff)
  value = time_diff + b / vk;
end

% c errors against time differences
time_differences = 2e-6:0.4e-7:8e-6;
errors_time_differences = zeros(size(time_differences)(2), 1);
time_differences_skews = zeros(size(time_differences)(2), 1);

for i = 1:size(errors_time_differences)(1)
  value = get_value_time_difference_known(b, vk, time_differences(i));

  time_differences_skews(i) = calculate_error_skew(error, value);
  minus_error = calculate_negative_error(error, value);
  plus_error = calculate_positive_error(error, value);
  % Multiply with b to rescale everything
  errors_time_differences(i) = (plus_error - minus_error) * b;
endfor

% c errors against weight percentage
weight_percentages = 0:0.1:20;
errors_weight_percentages = zeros(size(weight_percentages)(2), 1);
weight_percentages_skews = zeros(size(weight_percentages)(2), 1);

for i = 1:size(errors_weight_percentages)(1)
  value = get_value_v_known(b, calculate_speed_of_sound(set_temperature, weight_percentages(i)));

  weight_percentages_skews(i) = calculate_error_skew(error, value);
  minus_error = calculate_negative_error(error, value);
  plus_error = calculate_positive_error(error, value);
  % Multiply with b to rescale everything
  errors_weight_percentages(i) = (plus_error - minus_error) * b;
endfor

% c errors against rounding errors
errors = 0:1e-10:7e-8;
errors_from_rounding_errors = zeros(size(errors)(2), 1);
rounding_errors_skews = zeros(size(errors)(2), 1);

for i = 1:size(errors_from_rounding_errors)(1)
  value = get_value_v_known(b, calculate_speed_of_sound(set_temperature, set_weight_percentage));

  rounding_errors_skews(i) = calculate_error_skew(errors(i), value);
  minus_error = calculate_negative_error(errors(i), value);
  plus_error = calculate_positive_error(errors(i), value);
  % Multiply with b to rescale everything
  errors_from_rounding_errors(i) = (plus_error - minus_error) * b;
endfor


figure(1);
plot(time_differences, errors_time_differences, "linewidth", 8);
set(gca, "linewidth", 6, "fontsize", 22);
title(sprintf("c error due to rounding error of %.3e µs against time difference", error));
xlabel("time difference in s");
ylabel("delta c in m/s");

figure(2);
plot(time_differences, time_differences_skews, "linewidth", 8);
set(gca, "linewidth", 6, "fontsize", 22);
title(sprintf("c error skew (see notes) due to rounding error of %.3e µs against time difference", error));
xlabel("time difference in s");
ylabel("error skew");

figure(3);
plot(weight_percentages, errors_weight_percentages, "linewidth", 8);
set(gca, "linewidth", 6, "fontsize", 22);
title(sprintf("c error due to rounding error of %.3e µs against wp at T=%.2f °C", error, set_temperature));
xlabel("weight percentage in %");
ylabel("delta c in m/s");

figure(4);
plot(weight_percentages, weight_percentages_skews, "linewidth", 8);
set(gca, "linewidth", 6, "fontsize", 22);
title(sprintf("c error skew (see notes) due to rounding error of %.3e µs against time difference at T=%.2f °C", error, set_temperature));
xlabel("weight percentage in %");
ylabel("error skew");


figure(5);
plot(errors, errors_from_rounding_errors, "linewidth", 8);
set(gca, "linewidth", 6, "fontsize", 22);
title(sprintf("c error against a rounding error at wp=%.2f %% and T=%.2f °C", set_weight_percentage, set_temperature));
xlabel("rounding error of t in +-t");
ylabel("delta c in m/s");

figure(6);
plot(errors, rounding_errors_skews, "linewidth", 8);
set(gca, "linewidth", 6, "fontsize", 22);
title(sprintf("c error skew against a rounding error at wp=%.2f %% and T=%.2f °C", set_weight_percentage, set_temperature));
xlabel("rounding error of t in +-t");
ylabel("error skew");

