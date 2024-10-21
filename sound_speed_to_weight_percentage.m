% Obtained from https://iopscience.iop.org/article/10.1088/1742-6596/279/1/012023/pdf
% graphics_toolkit("fltk");

calculate_functions

function relative_error = calculate_relative_error(assumed, real)
  relative_error = (assumed / real) - 1;
endfunction

function [epsilon_sos, epsilon_wp, error_wp] = test_weight_percentage_backwards_calculation(target_temperature, target_wp)
  calculated_speed_of_sound = calculate_speed_of_sound(target_temperature, target_wp);
  [calculated_weight_percentage, epsilon_sos] = calculate_weight_percentage(target_temperature, calculated_speed_of_sound);

  epsilon_wp = abs(calculated_weight_percentage - target_wp);
  error_wp = calculate_relative_error(calculated_weight_percentage, target_wp);
endfunction

function conduct_random_tests(n)
  for i = 1:n
    target_temperature = rand() * 50;
    target_speed_of_sound = rand() * 25;
    [epsilon_sos, epsilon_wp, error_wp] = test_weight_percentage_backwards_calculation(target_temperature, target_speed_of_sound);
    printf("Attempt %d Error=%d T=%d WP=%d%s", i, error_wp,  target_temperature, target_speed_of_sound, "\n");
  endfor
endfunction

temperature = 20;
c_values_x = 1500:0.1:1600;
c_values_y = zeros(1, length(c_values_x));

for i = 1:length(c_values_x)
  c_values_y(i) = calculate_weight_percentage(temperature, c_values_x(i));
  printf("%.2f %% done\n", (i - 1) / length(c_values_x) * 100);
endfor

figure(1);
plot(c_values_x, c_values_y);
title(sprintf("Weight percentage against speed of sound c at T=%.2f", temperature));





