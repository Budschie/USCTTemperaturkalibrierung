% Evaluates the experiment using the CSV files
graphics_toolkit("fltk");

pkg load signal
pkg load image

calculate_functions

global GLOBAL_ORDER;
GLOBAL_ORDER = 50;

function plot_predicted_c_error(temperatures, time, measured_velocities, first_weight_percentage)
  predicted_sos = calculate_bulk_speed_of_sound(temperatures, repmat(first_weight_percentage, length(temperatures), 1));

  figure(4);
  hold on;
  set(gca, "linewidth", 6, "fontsize", 22);
  scatter(time, predicted_sos - measured_velocities);
  plot(time, predicted_sos - measured_velocities, "linewidth", 4);
  xlabel("t in min");
  ylabel("\Delta c in m/s (predicted - measured)");
  title("Difference due to different starting wp and change to said wp.");
endfunction

function evaluate_experiment_1
  sheet = csvread("EthanolEvaporationTest.csv");

  t_ende = sheet(3:end,13);
  T_oben = sheet(3:end,4);
  T_unten = sheet(3:end,5);
  v = sheet(3:end,8);

  [wp_oben, epsilon] = calculate_bulk_weight_percentage(T_oben, v);
  [wp_unten, epsilon] = calculate_bulk_weight_percentage(T_unten, v);

  figure(3);
  hold on;
  set(gca, "linewidth", 6, "fontsize", 22);
  scatter(t_ende, wp_oben);
  scatter(t_ende, wp_unten);
  plot(t_ende, wp_oben, "linewidth", 4);
  plot(t_ende, wp_unten, "linewidth", 4);
  ylabel("WP in %");
  xlabel("t in min");
  title("Weight percentage change across time");
  legend("WP mit T_o_b_e_n", "WP mit T_u_n_t_e_n");
  hold off;

  plot_predicted_c_error(T_oben, t_ende, v, wp_oben(2));
endfunction

% Signal is assumed to be in the last listed channel
function [time, signal, resolution, delay] = read_signal_csv(name, measurement_directory)
  path = sprintf("../%s/%s", measurement_directory, name);
  signal_output_file = csvread(path);

  signal_offset = size(signal_output_file, 2) - 5;

  delay = signal_output_file(2, 3 + signal_offset);
  resolution = signal_output_file(2, 4 + signal_offset);

  time_unsampled = signal_output_file(3:end, 1) * resolution + delay;
  signal_unsampled = signal_output_file(3:end, 2 + signal_offset);

  time = time_unsampled(1:2:end);
  signal = signal_unsampled(1:2:end);

  resolution *= 2;
endfunction

function [lowpassed_signal, lowpassed_time] = perform_lowpass(time, signal, sample_rate)
  global GLOBAL_ORDER;
  order = GLOBAL_ORDER * 2;
  % Delay always 10, WHY THOUGH?
  delay = order / 2;
  nominator = fir1(order, 5e6/sample_rate,'low');
  lowpassed_signal = filter(nominator, 1, signal)(1+delay:end);
  lowpassed_time = time(1:end-delay);
  % nfft = length(signal);
  % f = (-nfft/2:nfft/2-1) * sample_rate / nfft;
  % lowpassed_signal = butter(signal, lowpass_frequency);
endfunction

function [new_time, new_signal] = perform_interpft(time, signal, resampling_factor)
  new_length = length(time) * resampling_factor;

  new_time = interp1([1, length(signal) * resampling_factor], [time(1), time(end)], 1:1:length(signal) * resampling_factor);
  new_signal = interpft(signal, new_length);
endfunction

function max_t = acquire_max_t(name, measurement_directory)
  [time, signal, resolution, delay] = read_signal_csv(name, measurement_directory);
  [lowpassed_signal, lowpassed_time] = perform_lowpass(time, signal, 1 / resolution);
  [lowpassed_interpft_time, lowpassed_interpft_sig] = perform_interpft(lowpassed_time, lowpassed_signal, 40);
  [new_time, new_signal] = perform_interpft(time, signal, 20);

  [lowpassed_interpft_y, lowpassed_interpft_index] = max(abs(lowpassed_interpft_sig));
  lowpassed_x = lowpassed_interpft_time(lowpassed_interpft_index);

  [raw_interpft_y, raw_interpft_index] = max(abs(new_signal));
  raw_x = new_time(raw_interpft_index);

  % INTERESTING THING: Using the lowpassed interpft increases wp calculation accuracy in a range of e-13 or even x2 (sometimes makes it worse tho)

  max_t = raw_x;

  figure(10);
  hold on;
  plot(new_time, new_signal);
endfunction

function kelvin = to_kelvin(temp)
  kelvin = temp + 273.15;
endfunction

% Fit not correct/BS
function [new_v_k, new_b] = calc_vk(temp)
  b = 15.4790e-3;
  v_kc = 2552;

  k_1 = b / nthroot(20, 3);
  k_2 = v_kc / nthroot(20, 2);

  new_b = nthroot(temp, 3) * k_1;
  new_v_k = nthroot(temp, 2) * k_2;
endfunction

function evaluate_experiment_3(empty_prefix, full_prefix, index1, index2, temp_1, temp_2)
  measurement_directory = "Messergebnisse"

  empty_name_1 = sprintf("%s%d.csv", empty_prefix, index1);
  empty_name_2 = sprintf("%s%d.csv", empty_prefix, index2);
  full_name_1 = sprintf("%s%d.csv", full_prefix, index1);
  full_name_2 = sprintf("%s%d.csv", full_prefix, index2);
  % full_name = sprintf("%s%d.csv", full_prefix, index);

  max_t_empty_1 = acquire_max_t(empty_name_1, measurement_directory);
  max_t_empty_2 = acquire_max_t(empty_name_2, measurement_directory);
  max_t_full_1 = acquire_max_t(full_name_1, measurement_directory);
  max_t_full_2 = acquire_max_t(full_name_2, measurement_directory);


  width = 0.13;

  real_diff = (max_t_empty_1 - max_t_empty_2) / width;
  % real_diff = 1 / (width / max_t_empty_1 - width / max_t_empty_2);


  wp_real = 12.646006;

  expected_diff = (width / calculate_speed_of_sound(temp_1, wp_real) - width / calculate_speed_of_sound(temp_2, wp_real)) / width;

  printf("Expected diff of %.4e; real diff of %.4e\n", expected_diff, real_diff);

  b = 15.4790e-3;
  v_k = 2552;

  printf("Estimated c %.4f %.4f\n", calculate_speed_of_sound_relative(max_t_empty_1, max_t_full_1, b, v_k), calculate_speed_of_sound_relative(max_t_empty_2, max_t_full_2, b, v_k));

  figure(5);
  wps = [0:0.05:20];
  hold on;

  all_diffs = (width ./ calculate_bulk_speed_of_sound(repmat(temp_1, length(wps)), wps) - width ./ calculate_bulk_speed_of_sound(repmat(temp_2, length(wps)), wps)) ./ width;
  plot(wps, all_diffs);
  plot(wps, repmat(expected_diff, size(wps)));
  plot(wps, repmat(real_diff, size(wps)));
  legend("All diffs", "Expected diff", "Real diff");

  hold off;
endfunction


function wp = evaluate_experiment_2(index, temperature, measurement_directory, empty_prefix, full_prefix)
  b = 15.4790e-3;
  v_k = 2552;

  % [b, v_k] = calc_vk(temperature);

  empty_name = sprintf("%s%d.csv", empty_prefix, index);
  full_name = sprintf("%s%d.csv", full_prefix, index);

  max_t_empty = acquire_max_t(empty_name, measurement_directory);
  max_t_full = acquire_max_t(full_name, measurement_directory);

  v = calculate_speed_of_sound_relative(max_t_empty, max_t_full, b, v_k);
  [wp, sos_difference] = calculate_weight_percentage(temperature, v);

  printf("Index %d with T=%.2e°C and c=%.4e has estimated weight percentage of wp=%.3f%%, being off by %.3e m/s\n", index, temperature, v, wp, sos_difference);
endfunction

% printf(" --- Ethanol-Water evaporation experiment --- \n");

% wp1 = evaluate_experiment_2(1, 22.84, "Messung26.4", "filer", "fivol");
% wp2 = evaluate_experiment_2(2, 33.72, "Messung26.4", "filer", "fivol");
% wp_diff = wp2 - wp1;

% sample_wp_1 = calculate_speed_of_sound(30, wp1);
% sample_wp_2 = calculate_speed_of_sound(30, wp2);
% sample_wp_diff = sample_wp_2 - sample_wp_1;

% printf("WP diff of %.6f creates m/s unaccuracy of %.6e\n", wp_diff, sample_wp_diff);

% printf("\n\n --- Water test experiment --- \n");

% temps = [23.16, 23.25, 23.36];

% for i = 1:3
%   evaluate_experiment_2(i, temps(i), "WaterMeasuremen", "wal", "wav");
% endfor

function show_spectrogram_graph(signals, resolutions)
  figure(5);
  hold on;
  set(gca, "linewidth", 6, "fontsize", 22);

  for i = 1:size(signals)(2)
    padded_signal = padarray(signals(1:end, i), length(signals(1:end, i)) * 8);

    frequency = 1 / resolutions(i);
    ffted_signal = fft(padded_signal);
    fft_shifted = fftshift(ffted_signal);

    n = length(fft_shifted);

    frequency_x = (-n/2:n/2-1) * (frequency / n);
    power = (1 / (frequency * n)) * abs(fft_shifted).^2;

    plot(frequency_x, power, "linewidth", 1.7);
  endfor
endfunction

function show_graph_experiment_2
  global GLOBAL_ORDER;

  [time, signal, resolution, delay] = read_signal_csv("filer1.csv", "Messung26.4");
  [lowpassed_signal, lowpassed_time] = perform_lowpass(time, signal, 1 / resolution);
  [new_time, new_signal] = perform_interpft(time, signal, 20);
  [lowpassed_interpft_time, lowpassed_interpft_sig] = perform_interpft(lowpassed_time, lowpassed_signal, 20);

  % show_spectrogram_graph([signal(1:end-GLOBAL_ORDER), lowpassed_signal], [resolution, resolution]);
  show_spectrogram_graph([signal(1:end)], [resolution]);
  legend("Normal", "Lowpass");

  figure(6);
  hold on;
  set(gca, "linewidth", 6, "fontsize", 22);
  plot(time, signal, "linewidth", 4);
  plot(new_time, new_signal, "linewidth", 4);
  % plot(lowpassed_time, lowpassed_signal, "linewidth", 4);
  % plot(lowpassed_interpft_time, lowpassed_interpft_sig, "linewidth", 4);
  xlabel("t in s");
  ylabel("U in V");
  title("Oscilloscope signal of transducer");
  legend("Raw", "interpft (raw)", "lowpass", "interpft (lowpass)");
  hold off;

  [raw_interpft_y, raw_interpft_index] = max(new_signal);
  [lowpassed_interpft_y, lowpassed_interpft_index] = max(lowpassed_interpft_sig);

  raw_x = new_time(raw_interpft_index);
  lowpassed_x = lowpassed_interpft_time(lowpassed_interpft_index);

  printf("Raw interpft %.5e\nLowpassed interpft %.5e\nDelta %.5e\n", raw_x, lowpassed_x, raw_x - lowpassed_x);
endfunction

% show_graph_experiment_2
evaluate_experiment_3("al", "av", 1, 2, 33.86, 38.88);
% evaluate_experiment_3("l", "v", 2, 3, 28.22, 28.5);
