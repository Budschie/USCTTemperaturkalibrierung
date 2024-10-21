global DI_INDEX KI_INDEX LI_INDEX speed_of_sound_matrix;

DI_INDEX = 2;
KI_INDEX = 3;
LI_INDEX = 4;

speed_of_sound_matrix = dlmread("UltrasonicSpeedList.csv", ",");

function speed_of_sound = calculate_speed_of_sound_old(temperature, ethanol_weight_percentage)
  global DI_INDEX KI_INDEX LI_INDEX speed_of_sound_matrix;
  e_power_k_i = ethanol_weight_percentage .^ speed_of_sound_matrix(:,KI_INDEX);
  t_power_l_i = temperature .^ speed_of_sound_matrix(:,LI_INDEX);

  catted = cat(2, e_power_k_i, t_power_l_i, speed_of_sound_matrix(:,DI_INDEX));

  speed_of_sound_elements = prod(catted, 2);
  speed_of_sound = sum(speed_of_sound_elements);
end

function [time, voltage] = load_csv_file(path)
  read_file = csvread(path);

  start_time = read_file(2,3);
  increment_time = read_file(2,4);

  time = start_time + increment_time * read_file(3:2:end, 1);
  voltage = read_file(3:2:end, 2);
endfunction

% Must be evenly sampled
function [interpft_time, interpft_voltage] = perform_interpft(time, voltage, n)
  interpft_voltage = interpft(voltage, length(voltage) * n);
  interpft_time = (0:(length(voltage) * n -1)) / (length(voltage) * n - 1) * (time(end) - time(1)) + time(1);
endfunction

% What if dimensions are different? Dimension of wp will be chosen
function speed_of_sound = calculate_speed_of_sound(temperatures, wps)
  global DI_INDEX KI_INDEX LI_INDEX speed_of_sound_matrix;

  size_temps = size(temperatures);
  size_wps = size(wps);

  if size_temps != size_wps && size_temps != size(wps')
    error("Sizes of temp and wp must be the same");
  endif

  flattened_temps = vec(temperatures)';
  flattened_wps = vec(wps)';

  resulting = zeros(1, length(flattened_temps));

  for i=1:size(speed_of_sound_matrix)(1)
    wp_exp = flattened_wps .^ speed_of_sound_matrix(i, KI_INDEX);
    t_exp = flattened_temps .^ speed_of_sound_matrix(i, LI_INDEX);
    d_factor = speed_of_sound_matrix(i, DI_INDEX);
    resulting += wp_exp .* t_exp .* d_factor;
  endfor

  % wp_exp = flattened_wps .^ speed_of_sound_matrix(:, KI_INDEX);
  % t_exp = flattened_temps .^ speed_of_sound_matrix(:, LI_INDEX);
  % d_factor = speed_of_sound_matrix(:, DI_INDEX);

  % resulting = wp_exp .* t_exp .* d_factor;

  % speed_of_sound = reshape(sum(resulting, 1), size_wps);
  speed_of_sound = reshape(resulting, size_wps);
endfunction

function perform_test()
  a = [rand() * 10 + 10, rand() * 10 + 10; rand() * 10 + 10, rand() * 10 + 10];
  b = [rand() * 10 + 10, rand() * 10 + 10; rand() * 10 + 10, rand() * 10 + 10];

  resulting_new = calculate_speed_of_sound(a, b);

  for i = 1:2
    for j = 1:2
      if resulting_new(i,j) != calculate_speed_of_sound_old(a(i,j), b(i,j))
        error("Test failed");
      endif
    endfor
  endfor

  printf("Test successful\n");
endfunction

% calculate_speed_of_sound([10, 20; 10, 20], [10, 20; 10, 20])

function [resulting_weight_percentage, epsilon_sos] = calculate_weight_percentage(target_temperature, target_speed_of_sound)
  % Speed of sound strictly monotically increasing in +weight percentage direction
  % 0 <= wp <= 25
  current_lower_wp = 0;
  current_higher_wp = 25;
  max_iter = 1000;

  current_wp = 0;
  current_speed_of_sound = 0;

  resulting_weight_percentage = 0;
  for i = 1:max_iter
    current_wp = mean([current_lower_wp, current_higher_wp]);
    current_speed_of_sound = calculate_speed_of_sound(target_temperature, current_wp);

    % printf("Attempt %d w/ current sos of %d with wp of %d", i, current_speed_of_sound, current_wp);
    % Found value smaller than expected => discard anything below current_wp
    if current_speed_of_sound < target_speed_of_sound
      current_lower_wp = current_wp;
    endif

    if current_speed_of_sound > target_speed_of_sound
      current_higher_wp = current_wp;
    endif
  endfor

  resulting_weight_percentage = current_wp;
  epsilon_sos = abs(target_speed_of_sound - current_speed_of_sound);
endfunction


function speed_of_sounds = calculate_bulk_speed_of_sound(temperatures, weight_percentages)
  speed_of_sounds = calculate_speed_of_sound(temperatures, weight_percentages);
end

function [weight_percentages, epsilon_sos] = calculate_bulk_weight_percentage(temperatures, speed_of_sounds)
  weight_percentages = zeros(length(temperatures), 1);
  epsilon_sos = zeros(length(temperatures), 1);

  for i = 1:length(temperatures)
    [weight_percentages(i), epsilon_sos(i)] = calculate_weight_percentage(temperatures(i), speed_of_sounds(i));
  endfor
end

function speed_of_sound = calculate_speed_of_sound_relative(t_empty, t_full, b, v_k)
  speed_of_sound = b ./ (t_empty - t_full + b ./ v_k);
endfunction

function speed_of_sound = calculate_reverse_speed_of_sound_relative(t_empty, t_full, b, v_k)
  speed_of_sound = b ./ (t_full - t_empty + b ./ v_k);
endfunction


