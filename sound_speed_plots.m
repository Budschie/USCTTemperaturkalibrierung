% Obtained from https://iopscience.iop.org/article/10.1088/1742-6596/279/1/012023/pdf
graphics_toolkit("qt");

% Note: This file has been superceeded by c_plot_2.m

calculate_functions

% Parameters for the sample area
% temperature_min = 18;
% temperature_max = 22.0;
% temperature_min_max_delta = temperature_max - temperature_min;

% alcohol_weight_percentage_min = 0.0;
% alcohol_weight_percentage_max = 20.0;
% alcohol_weight_percentage_min_max_delta = alcohol_weight_percentage_max - alcohol_weight_percentage_min;

% resolution_temperature = (temperature_min_max_delta) * 2;
% resolution_alcohol_weight_percentage = (alcohol_weight_percentage_min_max_delta) * 2;

% X: Temperature; Y: Weight percentage; Z: Speed of sound
%% evaluated_values = zeros(resolution_temperature + 1, resolution_alcohol_weight_percentage + 1);

% Create a grid of samples of the function calculate_speed_of_sound
%% for x = 0:resolution_temperature
%%  current_temperature = (x / resolution_temperature) * temperature_min_max_delta + temperature_min;

%%  for y = 0:resolution_alcohol_weight_percentage
%%    current_alcohol_weight_percentage = (y / resolution_alcohol_weight_percentage) * alcohol_weight_percentage_min_max_delta + alcohol_weight_percentage_min;

    % Sample the function and store its result
%%    evaluated_values(x + 1, y + 1) = calculate_speed_of_sound(current_temperature, current_alcohol_weight_percentage);
%%  endfor
%% endfor

% Generate mesh grid

expected_weight_perc = 13.0;

global min_temp max_temp;

min_temp = 25;
max_temp = 40;
temp_resolution = 0.05;

min_weight_perc = 0;
max_weight_perc = 20;
weight_perc_resolution = 0.025;

% Calculate the deviation from an expected speed of sound value given the weight_percentage and a temperature due to inaccuracies when measuring the temperature

% USED IN A WRONG WAY!!! USE A (CURRENTLY NOT IMPLEMENTED) FUNCTION THAT TAKES IN A TEMPERATURE UNCERTAINTY VALUE INSTEAD!!!
% (ALTHOUGH, MEDIUM IS NOT UNIFORM SO IT MIGHT BE LEGIT WHAT I AM DOING HERE)
function temperature_caused_sound_deviation = calculate_temperature_caused_sound_deviation(weight_percentage)
  global min_temp max_temp;
  expected_temperature = (min_temp + max_temp) * 0.5;

  actual_sos = calculate_speed_of_sound(expected_temperature, weight_percentage);
  temperature_caused_sound_deviation = max(abs(actual_sos - calculate_speed_of_sound(min_temp, weight_percentage)),
    abs(actual_sos - calculate_speed_of_sound(max_temp, weight_percentage)));
end

##function weight_percentage_caused_sound_deviation_helper = calculate_weight_percentage_caused_sound_deviation_helper(expected_wp, expected_speed_of_sound, weight_perc_deviation)
##    global min_temp max_temp;
##
##    sos_min_temp = calculate_speed_of_sound(min_temp, expected_wp + weight_perc_deviation);
##    sos_max_temp = calculate_speed_of_sound(max_temp, expected_wp + weight_perc_deviation);
##
##    weight_percentage_caused_sound_deviation_helper = max(abs(expected_speed_of_sound - sos_min_temp), abs(expected_speed_of_sound - sos_max_temp));
##end
##
##% Calculate the deviation from an expected speed of sound at a specified weight percentage and temperature due to changing temperature and changing weight percentage
##function weight_percentage_caused_sound_deviation = calculate_weight_percentage_caused_sound_deviation(weight_percentage, weight_percentage_deviation)
##  global min_temp max_temp;
##  expected_temperature = (min_temp + max_temp) * 0.5;
##  expected_value = calculate_speed_of_sound(expected_temperature, weight_percentage);
##
##  % TODO: FIX THIS FUNCTION
##
##  sos_deviation_min_wp = calculate_weight_percentage_caused_sound_deviation_helper(weight_percentage, expected_value, -weight_percentage_deviation);
##  sos_deviation_max_wp = calculate_weight_percentage_caused_sound_deviation_helper(weight_percentage, expected_value, +weight_percentage_deviation);
##
##  weight_percentage_caused_sound_deviation = max(sos_deviation_max_wp, sos_deviation_min_wp);
##  % weight_percentage_caused_sound_deviation = max(abs(calculate_temperature_caused_sound_deviation(weight_percentage + weight_percentage_deviation)),
##  %  abs(calculate_temperature_caused_sound_deviation(weight_percentage - weight_percentage_deviation)));
##end

function error_helper = calculate_speed_of_sound_deviation_from_w_err_w_off_helper(temperature, target_weight_percentage, w_err, w_off)
  global min_temp max_temp;
  expected_temp = (min_temp + max_temp) * 0.5;
  real_world_sos = calculate_speed_of_sound(temperature, target_weight_percentage + w_off);
  projected_sos = calculate_speed_of_sound(expected_temp, target_weight_percentage + w_off + w_err);

  error_helper = abs(real_world_sos - projected_sos);
end

% Where's error incorporated by temperature stuff?
function error = calculate_speed_of_sound_deviation_from_w_err_w_off(target_weight_percentage, w_err, w_off)
  global min_temp max_temp;
  min_sos = calculate_speed_of_sound_deviation_from_w_err_w_off_helper(min_temp, target_weight_percentage, w_err, w_off);
  max_sos = calculate_speed_of_sound_deviation_from_w_err_w_off_helper(max_temp, target_weight_percentage, w_err, w_off);
  % real_world_sos_deviation = calculate_temperature_caused_sound_deviation(target_weight_percentage + w_off);
  % projected_sos_deviation = calculate_temperature_caused_sound_deviation(target_weight_percentage + w_off + w_err);
  error = max(min_sos, max_sos);
end

[X, Y] = meshgrid(min_temp:temp_resolution:max_temp, min_weight_perc:weight_perc_resolution:max_weight_perc);
[X2, Y2] = meshgrid(min_temp:temp_resolution:max_temp, (min_weight_perc - expected_weight_perc):weight_perc_resolution:(max_weight_perc - expected_weight_perc));
[X3, Y3] = meshgrid(-5:0.2:5, -5:0.2:5);

x_size3 = size(X3)(1);
y_size3 = size(X3)(2);
evaluated_diff_values = zeros(x_size3, y_size3);
baseline_diff_values = zeros(x_size3, y_size3);

baseline_diff_value = calculate_temperature_caused_sound_deviation(0);

for x2 = 1:x_size3
  for y2 = 1:y_size3
    evaluated_diff_values(x2, y2) = calculate_speed_of_sound_deviation_from_w_err_w_off(expected_weight_perc, X3(x2, y2), Y3(x2, y2));
    baseline_diff_values(x2, y2) = baseline_diff_value;
  endfor
endfor

function graph_weight_percentage_at_temperature(wp_range, temperature)
  wp_size = (size(wp_range))(2)(1);
  results = zeros(wp_size, 1);

  for x = 1:wp_size
    results(x) = calculate_speed_of_sound(temperature, wp_range(x)(1));
  endfor

  plot(wp_range, results);
endfunction

% Ethanol deviation X, Ethanol
##EDX = 0:0.01:20;
##
##edy_size = size(EDX)(2);
##edy_created = zeros(edy_size)(1);
##edy_created_2 = zeros(edy_size)(1);
##
##ref_value = calculate_temperature_caused_sound_deviation(0);
##
##for x = 1:edy_size
##  edy_created(x) = calculate_weight_percentage_caused_sound_deviation(20, EDX(x));
##  edy_created_2(x) = ref_value;
##end

x_size = size(X)(1);
y_size = size(X)(2);
evaluated_values = zeros(x_size, y_size);
evaluated_values_deviations = zeros(x_size, y_size);
error_created = zeros(x_size, y_size);

for x = 1:x_size
  for y = 1:y_size
    expected_speed_of_sound = calculate_speed_of_sound(X2(x,y), expected_weight_perc);
    actual_speed_of_sound = calculate_speed_of_sound(X2(x, y), expected_weight_perc + Y2(x,y));
    error_created(x, y) = abs(actual_speed_of_sound - expected_speed_of_sound);

    % Calculate the speed of sound at a temperature of 20 degrees for the given alcohol density
    target_speed_of_sound = calculate_speed_of_sound((min_temp + max_temp) / 2.0, Y(x,y));
    % Calculate the speed of sound at the given temperature and alcohol density
    evaluated_speed_of_sound = calculate_speed_of_sound(X(x,y), Y(x,y));
    evaluated_values(x, y) = evaluated_speed_of_sound;
    % Absolute error of the calculated speed of sound and the target speed of sound
    evaluated_values_deviations(x, y) = abs(evaluated_speed_of_sound - target_speed_of_sound);
  endfor
endfor

changing_temperature = min_temp:temp_resolution:max_temp;
sos_for_temp = zeros(size(changing_temperature)(2));

for i = 1:size(sos_for_temp)(1)
  sos_for_temp(i) = calculate_speed_of_sound(changing_temperature(i), expected_weight_perc);
endfor

figure(11);
plot(changing_temperature, sos_for_temp);
title("Speed of sound against temperature for a given weight percentage");
xlabel("T in °C");
ylabel("c in m/s");

figure(1);
pcolor(X, Y, evaluated_values);
title("Alcohol weight percentage and temperature to speed of sound");
xlabel("Temperature in °C");
ylabel("Alcohol weight percentage in percent");
zlabel("Speed of sound in m/s");
colorbar();

figure(2);
pcolor(X, Y, evaluated_values_deviations);
title("Deviations of Speed of sound from temperature measured at 20°C with the given alcohol weight percentage");
xlabel("Temperature in °C");
ylabel("Alcohol weight percentage in percent");
zlabel("Speed of sound deviation in m/s");
colorbar();


figure(3);
mean_value_for_alcohol_density = mean(evaluated_values_deviations, 2);
plot(Y, mean_value_for_alcohol_density);
title("Alcohol weight percentage against mean speed of sound deviation from a sampled speed of sound at 20°C");
xlabel("Alcohol weight percentage");
ylabel("Average deviation from speed of sound at 20°C in m/s");

% Graph difference between expected sound speed and actual sound speed
% X axis: temperature
% Y axis: deviation of volume
% Z axis: sound speed difference
% figure(4);
% Goal: Determine how the accuracy would change according to a deviation from the desired mixing ratio
figure(4);
pcolor(X2, Y2, error_created);
title("Difference between predicted speed of sound and actual speed of sound due to uncertainty in the ethanol water mixing process");
xlabel("Temperature in °C");
ylabel("Alcohol weight percentage difference");
zlabel("Difference of speed of sound in m/s");
colorbar();

figure(5);
mean_value_for_prediction_vs_reality = mean(error_created, 2);
plot(Y2, mean_value_for_prediction_vs_reality);
title("Mean difference between predicted speed of sound and actual speed of sound");
xlabel("Alcohol weight percentage difference");
ylabel("Mean difference of speed of sound in m/s");

figure(6);
% Calculate difference between the speed of sound at the smallest and highest temperature
% Actually weight percentage
tempstuff = min_weight_perc:weight_perc_resolution:max_weight_perc;
minmaxdiff = zeros(size(tempstuff)(2));

for i = 1:size(minmaxdiff)
  minmaxdiff(i) = calculate_speed_of_sound(min_temp, tempstuff(i)) - calculate_speed_of_sound(max_temp, tempstuff(i));
endfor

plot(tempstuff, minmaxdiff);
title("Difference between speed of sound at min_temp and max_temp using a given alcohol weight percentage");
xlabel("Weight percentage");
ylabel("Deviation of speed of sound");

[min_value, min_index] = min(abs(minmaxdiff));
printf("Min difference is %d, min weight percentage %.*d", minmaxdiff(min_index(1)), 6, tempstuff(min_index(1)));

figure(7);
% X-Axis: Err
% Y-Axis: W_Off

hold on

surface(X3, Y3, evaluated_diff_values);
title("Speed of sound errors due to wrong mixing and wrong calibration");
xlabel("Weight percentage offset due to wrong calibration");
ylabel("Weight percentage offset due to wrong mixing");

surface(X3, Y3, baseline_diff_values, "FaceColor", [1, 0, 0]);

hold off

figure(8);
graph_weight_percentage_at_temperature((expected_weight_perc - 1):0.5:(expected_weight_perc + 1), 20);
##plot(EDX, edy_created);
##title("Max speed of sound deviation occouring within given temperature range and for a given weight percentage deviation");
##xlabel("Weight percentage deviation in absolute percent");
##ylabel("Max speed of sound deviation in m/s");
##
##hold on
##plot(EDX, edy_created_2);
##hold off

"Expected"
