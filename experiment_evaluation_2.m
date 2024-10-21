% Evaluates the experiment using the CSV files
graphics_toolkit("fltk");

global old_fit_procedure;
old_fit_procedure=false;

pkg load signal;

global save_plots;

save_plots="doc";

if save_plots == "doc"
  set(0, "defaultlinelinewidth", 1);
elseif save_plots == "presentation"
  set(0, "defaultlinelinewidth", 3);
else
  set(0, "defaultlinelinewidth", 3);
endif

global documentation_path;
documentation_path = "/home/budschie/Dokumente/HectorCloud/HectorData/Abschlussdokumentation";

calculate_functions;

pkg load image;

global noise_uncertainty_floor;
global noise_precision_margin;

noise_uncertainty_floor = 0.0004;
noise_precision_margin = 0.0001;

% CDF and PDF of continous uniform distribution math functions (not the octave ones) from https://en.wikipedia.org/wiki/Continuous_uniform_distribution

% Smallest theoretically possible signal strength at this location
function h_min = get_h_min(index, values)
  global noise_uncertainty_floor;
  h_min = values(index) - noise_uncertainty_floor;
endfunction

% Biggest theoretically possible signal strength at this location
function h_max = get_h_max(index, values)
  global noise_uncertainty_floor;
  h_max = values(index) + noise_uncertainty_floor;
endfunction

% Probability density, that the height y is achieved at the point index
function height_pd = get_height_pd(y, index, values)
  global noise_uncertainty_floor;
  if_expression = abs(y - values(index)) <= noise_uncertainty_floor;
  one_over_expression = 1 / (get_h_max(index, values) - get_h_min(index, values));

  height_pd = transpose(if_expression * one_over_expression);
endfunction

% Cumulated probability, that the height achieved is less then or equal leq_y at index
function height_cdf_capped = get_height_cdf(leq_y, indices, values)
  nominator = transpose(leq_y) - get_h_min(indices, values);
  denominator = (get_h_max(indices, values) - get_h_min(indices, values));
  height_cdf_capped = max(0, min(1, nominator ./ denominator));
endfunction

% Finds the first and the last value of t, which could theoretically due to noise also be the location of the peak
% Only important for optimization
function [min_t_bound, max_t_bound] = get_t_bounds(values)
  global noise_uncertainty_floor;
  global noise_precision_margin;
  % Find biggest value
  biggest = max(values)(1);

  % Subtract the noise uncertainty floor, as the real sampled value
  % could theoretically be below the current value
  % Strictly speaking, this probably neutralises with another term
  biggest -= noise_uncertainty_floor;
  biggest -= noise_precision_margin;

  % Create absolute difference to all values
  % Apply a margin to avoid floating point precision errors
  diff = (values + noise_uncertainty_floor) - biggest;

  min_t_bound = length(diff) + 1;
  max_t_bound = 0;

  for i=1:length(diff)
    if diff(i) >= 0
      max_t_bound = max(max_t_bound, i);
      min_t_bound = min(min_t_bound, i);
    endif
  endfor
endfunction

function [time, voltage] = load_csv_file(path)
  read_file = csvread(path);

  start_time = read_file(2,3);
  increment_time = read_file(2,4);

  time = start_time + increment_time * read_file(3:2:end, 1);
  voltage = read_file(3:2:end, 2);
endfunction

function [time, voltage] = interpft_filter(input_time, input_voltage, resampling_factor)
  time = transpose(interp1(1:length(input_time), input_time, 1:1/resampling_factor:length(input_time)));
  voltage = real(interpft(input_voltage, length(time)));
endfunction

function p_peak_prob_result = get_peak_prob(y, index, values, min_index, max_index)
  p_peak_prob = get_height_pd(y, index, values);

  % for i=min_index:max_index
    % if i != index
      % p_peak_prob *= get_height_cdf(y, i, values);
    % endif
  % endfor

  all_heights = get_height_cdf(y, min_index:max_index, values);
  all_heights(index - min_index + 1) = 1;


  cumulative_product = cumprod(all_heights, 1)(end,:);
  % over_cdf = 1 ./ get_height_cdf(y, index, values);
  p_peak_prob_result = prod([p_peak_prob; cumulative_product], 1);
  % p_peak_prob_result = p_peak_prob * cumulative_product / get_height_cdf(y, index, values);
endfunction

function [all_peak_probs, heights] = get_all_peak_probs(samples, index, values, min_index, max_index)
  global noise_uncertainty_floor;
  heights = (((-samples:1:samples) / samples) * noise_uncertainty_floor) + values(index);

  all_peak_probs = get_peak_prob(transpose(heights), index, values, min_index, max_index);

  % all_peak_probs = zeros(length(heights), 1);

  % for i=1:length(all_peak_probs)
    % all_peak_probs(i) = get_peak_prob(heights(i), index, values, min_index, max_index);
  % endfor
  all_peak_probs = get_peak_prob(transpose(heights), index, values, min_index, max_index);
endfunction

function resulting_pdf = get_all_peak_probs_integrated(all_peak_probs, heights)
  dh = heights(2) - heights(1);

  resulting_pdf = sum(all_peak_probs);
  resulting_pdf *= dh;
endfunction

function [density_mu, density_sigma] = get_pdf_stats(pdf, x_values)
  % Since heights are element of real numbers, but time is quantized, a PDF is not present anymore; instead there is a discrete probability
  dh = x_values(2) - x_values(1);

  ref_val = sum(pdf);
  printf("PDF validity of %.4f\n", ref_val);


  density_mu = sum(prod([pdf; x_values]));

  density_sigma = 0;

  for i = 1:length(pdf)
    to_sq = (x_values(i) - density_mu);
    density_sigma += pdf(i) * to_sq * to_sq;
  endfor

  density_sigma = sqrt(density_sigma);
endfunction

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

function [time, voltage] = moving_average_filter(input_time, input_voltage, window)
  filtered_voltage_pre = filter(ones(window,1)/window, 1, input_voltage); %# moving average
  delay = (window - 1) / 2;

  time = input_time(1:end-delay);
  voltage = filtered_voltage_pre(1+delay:end);
endfunction

function [time, voltage] = lowpass_filter(input_time, input_voltage)
  nyquist_frequency = 1/(input_time(2) - input_time(1)) / 2;

  cutoff=0.33*nyquist_frequency;

  n=1;

  [filter_butter_coeff, filter_butter_den_coeff] = butter(n, cutoff/nyquist_frequency);
  % [delay, frequencies] = grpdelay(filter_butter_coeff, filter_butter_den_coeff, 512, 1/(input_time(2) - input_time(1)));
  [delay, samples] = grpdelay(filter_butter_coeff, filter_butter_den_coeff);

  avg_delay = round(mean(delay));

  time = input_time(1:end-avg_delay);
  voltage = filter(filter_butter_coeff, filter_butter_den_coeff, input_voltage)((1+avg_delay):end);
endfunction

% Displays the waveform of the nth waveform file of the water experiment data
function p_an(n)
  global noise_uncertainty_floor;
  global noise_precision_margin;

  [time, voltage] = load_csv_file(sprintf("WaterFit/%d.csv", n));

  figure(1000);
  plot(time, voltage);
  title(sprintf("Measurement %d", n));
  xlabel("$t$ in s");
  ylabel("$U$ in V");

  figure(1);
  clf;
  plot(time, voltage);
  hold on;
  [interpft_time_nt, interpft_voltage] = interpft_filter(time, voltage, 10);
  interpft_time = transpose(interpft_time_nt);
  plot(interpft_time, interpft_voltage);

  % [moving_average_time, moving_average_voltage] = moving_average_filter(time, voltage, 13);
  % plot(moving_average_time, moving_average_voltage);

  % [interpft_time_nt_ma, interpft_voltage_ma] = interpft_filter(moving_average_time, moving_average_voltage, 10);
  % interpft_time_ma = transpose(interpft_time_nt_ma);
  % plot(interpft_time_ma, interpft_voltage_ma);

  [butter_time, butter_voltage] = lowpass_filter(time, voltage);
  plot(butter_time, butter_voltage);

  [interpft_time_nt_bt, interpft_voltage_bt] = interpft_filter(butter_time, butter_voltage, 10);
  interpft_time_bt = transpose(interpft_time_nt_bt);
  plot(interpft_time_bt, interpft_voltage_bt);

  % figure(5);
  % clf;
  % show_spectrogram_graph(interpft_voltage, interpft_time(2) - interpft_time(1));
  % show_spectrogram_graph(interpft_voltage_bt, interpft_time_bt(2) - interpft_time_bt(1));

  hold off;
  title(sprintf("Measurement %d", n));
  legend("Normal", "Interpft", "Butter", "Butter + interpft"); % "Moving average", "Interpft + moving average",

  figure(21);
  clf;
  hold on;
  plot(time, voltage, "o", "linewidth", 0.25);
  plot(interpft_time, interpft_voltage);
  title(sprintf("Measurement %d", n));
  xlabel("$t$ in s");
  ylabel("$U$ in V");

  figure(22);
  clf;
  hold on;
  plot(butter_time, butter_voltage, "o", "linewidth", 0.25);
  plot(interpft_time_bt, interpft_voltage_bt);
  title(sprintf("Measurement %d", n));
  xlabel("$t$ in s");
  ylabel("$U$ in V");

  figure(2);
  plot(time, voltage + noise_uncertainty_floor);
  hold on;
  % Noise precision margin is a delta value to avoid floating point precision errors later on
  plot(time, repmat(max(voltage)(1) - noise_uncertainty_floor - noise_precision_margin, length(time), 1));
  hold off;

  used_voltage = transpose(interp1(time, voltage, interpft_time));
  used_time = transpose(interpft_time);

  [minb, maxb] = get_t_bounds(used_voltage);

  figure(3);
  [_, max_voltage_index] = max(used_voltage);

  all_integrated_values = zeros(maxb - minb + 1, 1);


  for i=minb:maxb
    [peak_probs, heights] = get_all_peak_probs(250, i, used_voltage, minb, maxb);

    all_integrated_values(i - minb + 1) = get_all_peak_probs_integrated(peak_probs, heights);
    % plot(heights, peak_probs);

    % input(sprintf("Currently viewn sample is %d\n", i));
  endfor

  plot(used_time(minb:maxb), all_integrated_values);
  [mu, sigma] = get_pdf_stats(all_integrated_values, used_time(minb:maxb));
  printf("µ=%.6e; sigma=%.6e\n", mu, sigma);
  xlabel("t_p_e_a_k in s");
  ylabel("P_p_e_a_k");
  title("Probability of peaks against time");

  printf("Earliest possible peak at t=%.4e, latest possible at t=%.4e\n", used_time(minb), used_time(maxb));

  % figure(5);
  % clf;
  % show_spectrogram_graph([voltage(1:end-6), moving_average_voltage], repmat(time(2) - time(1), 2, 1));
endfunction

% Analyses the waveform, interpolates it using interpft and then returns the peak probabilities and the time of said peaks
function [time, probabilities] = analyse_waveform(prefix, number)
  global noise_uncertainty_floor;
  global noise_precision_margin;

  % Loading CSV
  [input_time, voltage] = load_csv_file(sprintf("%s/%d.csv", prefix, number));

  % Interpft
  [interpft_time_nt, interpft_voltage] = interpft_filter(input_time, voltage, 10);
  interpft_time = transpose(interpft_time_nt);

  used_voltage = transpose(interp1(input_time, voltage, interpft_time));
  used_time = transpose(interpft_time);

  % Instantly rejecting datapoints with 0% probability of being a peak
  [minb, maxb] = get_t_bounds(used_voltage);

  all_integrated_values = zeros(maxb - minb + 1, 1);

  % For every data point...
  for i=minb:maxb
    % Retrieve the set of heights and their corresponding probability of being a peak
    [peak_probs, heights] = get_all_peak_probs(250, i, used_voltage, minb, maxb);
    % Integrate over height as it is continous to get the total probability of the datapoint being a peak
    all_integrated_values(i - minb + 1) = get_all_peak_probs_integrated(peak_probs, heights);
  endfor

  plot(used_time(minb:maxb), all_integrated_values);

  time = used_time(minb:maxb);
  probabilities = all_integrated_values;
endfunction

function [filtered_time, filtered_voltage] = perform_filtering(input_time, input_voltage)
  [butter_filtered_time, butter_filtered_voltage] = lowpass_filter(input_time, input_voltage);

  [interpft_time_nt, interpft_voltage] = interpft_filter(butter_filtered_time, butter_filtered_voltage, 10);

  filtered_time = interpft_time_nt';
  filtered_voltage = interpft_voltage;
endfunction

function peak_time = perform_peak_detection(time, voltage)
  [max_value, max_index] = max(voltage);
  peak_time = time(max_index);
endfunction

function peak_time = perform_filtering_peak_detection(prefix, number)
  [time, voltage] = load_csv_file(sprintf("%s/%d.csv", prefix, number));

  [filtered_time, filtered_voltage] = perform_filtering(time, voltage);
  peak_time = perform_peak_detection(filtered_time, filtered_voltage);
endfunction

function sos = get_speed_of_sound(prefix, empty_number, full_number, temperature)
  peak_empty = perform_filtering_peak_detection(prefix, empty_number);
  peak_full = perform_filtering_peak_detection(prefix, full_number);

  sos = calculate_reverse_speed_of_sound_relative(peak_empty, peak_full, 11.7408e-3, get_water_sos(temperature));
endfunction

% https://www.sciencedirect.com/science/article/abs/pii/S030156299800091X
function water_sos = get_water_sos(temperature_c)
  water_sos = 1404.3 + 4.7 * temperature_c - 0.04 * temperature_c * temperature_c;
endfunction

function expected_value = show_surface_plot_probabilities(prefix, empty_number, full_number, temperature)
  [time_empty, probabilities_empty] = analyse_waveform(prefix, empty_number);
  [time_full, probabilities_full] = analyse_waveform(prefix, full_number);

  [time_empty_grid, time_full_grid] = meshgrid(time_empty, time_full);
  [probabilities_empty_grid, probabilities_full_grid] = meshgrid(probabilities_empty, probabilities_full);

  x_rep = size(time_empty_grid)(1);
  y_rep = size(time_empty_grid)(2);

  % p_an(empty_number);

  speed_of_sounds = calculate_reverse_speed_of_sound_relative(time_empty_grid, time_full_grid, repmat(11.7408e-3, x_rep, y_rep), repmat(get_water_sos(temperature), x_rep, y_rep));

  expected_value = sum(sum(speed_of_sounds .* probabilities_empty_grid .* probabilities_full_grid));
  standard_deviation = sqrt(sum(sum(probabilities_empty_grid .* probabilities_full_grid .* (speed_of_sounds - expected_value) .* (speed_of_sounds - expected_value))));

  printf("Expected value of %.4e m/s; sigma %.4e\n", expected_value, standard_deviation);

  figure(60);
  clf;
  surface(time_empty_grid, time_full_grid, speed_of_sounds);
  title("Speed of sounds against peak times");
  colorbar;

  figure(61);
  clf;
  surf(time_empty_grid, time_full_grid, probabilities_empty_grid .* probabilities_full_grid);
  title("Combined peak probabilities");
  colorbar;
  xlabel("t_e_m_p_t_y in s");
  ylabel("t_f_u_l_l in s");
  zlabel("P of peak combination");
endfunction

% Reads in an index file of the experiment, loads the different waveforms
% and extracts the peak probs. The prefix is the directory containing the waveforms.
function create_water_fit(prefix, index_file)
  global old_fit_procedure;

  index_file = csvread(index_file);

  temperatures = zeros(length(index_file) / 2, 1);
  t_diff = zeros(length(index_file) / 2, 1);
  speed_of_sounds = zeros(length(index_file) / 2, 1);

  for i = 0:length(index_file)/2 - 1
    f_1 = index_file(i * 2 + 1, 1);
    f_2 = index_file(i * 2 + 2, 1);
    t_1 = index_file(i * 2 + 1, 3);
    t_2 = index_file(i * 2 + 2, 3);

    id_1 = index_file(i * 2 + 1, 2);
    id_2 = index_file(i * 2 + 2, 2);

    f_empty = "";
    f_full = "";
    t_empty = 0;
    t_full = 0;

    % 0 means empty, 1 means full

    if id_1 == 0 && id_2 == 1
      f_empty = f_1;
      f_full = f_2;
      t_empty = t_1;
      t_full = t_2;
    elseif id_2 == 0 && id_1 == 1
      f_empty = f_2;
      f_full = f_1;
      t_empty = t_2;
      t_full = t_1;
    else
      error("malformed data in index file");
    endif

    temperatures(i + 1) = mean([t_1, t_2]);

    if old_fit_procedure
      speed_of_sounds(i + 1) = show_surface_plot_probabilities(prefix, f_empty, f_full, mean([t_1, t_2]));
    else
      speed_of_sounds(i + 1) = get_speed_of_sound(prefix, f_empty, f_full, mean([t_1, t_2]));
    endif
    t_diff(i + 1) = t_2 - t_1;
  endfor

  figure(8);
  clf;
  plot(temperatures, speed_of_sounds);
  %hold on;
  % plot(temperatures, t_diff);
  % hold off;
  % ylim([0, max(speed_of_sounds) + 100]);

  title("EpoTek material speed of sound against temperature");
  xlabel("T in °C");
  ylabel("c in m/s");
  % matlab2tikz("epo_tek_c_vs_T.tex", "showInfo", false, 'height', '\fheight', 'width', '\fwidth');
endfunction

function create_signal_figure
  p_an(16);
  figure(1000);
  % print("/home/budschie/Dokumente/HectorCloud/HectorData/FinalerVortragII/Images/Signals/raw.svg", "-dsvg", "-S1920,1080");
  figure(1);
  % print("/home/budschie/Dokumente/HectorCloud/HectorData/FinalerVortragII/Images/Signals/multiple.svg", "-dsvg", "-S1920,1080");
endfunction

function create_peak_combination_figure
  show_surface_plot_probabilities("WaterFit", 21, 20, mean([37.22, 37.32]));
  figure(61);
  shading interp;

  path_to_print_to="/home/budschie/Dokumente/HectorCloud/HectorData/FinalerVortragII/Images/VerticalError/Real";

  print(sprintf("%s/3d_vertical_error.png", path_to_print_to), "-dpng", "-S1920,1080");

  p_an(20);
  figure(3);
  print(sprintf("%s/2d_full.png", path_to_print_to), "-dpng", "-S1920,1080");

  p_an(21);
  figure(3);
  print(sprintf("%s/2d_empty.png", path_to_print_to), "-dpng", "-S1920,1080");
endfunction

function create_water_fit_figures
  global save_plots;

  create_water_fit("WaterFit", "WaterFit/Index.csv");

  figure(8);

  if strcmp(save_plots, "doc")
    global documentation_path;
    matlab2tikz(sprintf("%s/parts/ergebnisse/signal_processing/ethanol_fit.tex", documentation_path), "showInfo", false, 'height', '\fheight', 'width', '\fwidth');
  elseif strcmp(save_plots, "presentation")
    print("/home/budschie/Dokumente/HectorCloud/HectorData/FinalerVortragII/Images/VerticalError/Real/epo_tek_fit.png", "-dpng", "-S1920,1080");
  endif
endfunction

function display_single_signal()
  global save_plots;

  p_an(16);
  figure(1000);

  x_zoom_range=[4.1e-05, 4.14e-05];
  y_zoom_range=[0.0105, 0.015];

  if strcmp(save_plots, "doc")
    global documentation_path;
    matlab2tikz(sprintf("%s/parts/material_und_methoden/signal_processing/signal.tex", documentation_path), "showInfo", false, 'height', '\fheight', 'width', '\fwidth');
    % axis(y_zoom_range, x_zoom_range);
    xlim(x_zoom_range);
    ylim(y_zoom_range);
    matlab2tikz(sprintf("%s/parts/material_und_methoden/signal_processing/signal_zoomed.tex", documentation_path), "showInfo", false, 'height', '\fheight', 'width', '\fwidth');

    figure(21);
    matlab2tikz(sprintf("%s/parts/ergebnisse/signal_processing/signal_raw.tex", documentation_path), "showInfo", false, 'height', '\fheight', 'width', '\fwidth');
    xlim(x_zoom_range);
    ylim(y_zoom_range);
    matlab2tikz(sprintf("%s/parts/ergebnisse/signal_processing/signal_raw_zoomed.tex", documentation_path), "showInfo", false, 'height', '\fheight', 'width', '\fwidth');

    figure(22);
    matlab2tikz(sprintf("%s/parts/ergebnisse/signal_processing/signal_lp.tex", documentation_path), "showInfo", false, 'height', '\fheight', 'width', '\fwidth');
    xlim(x_zoom_range);
    ylim(y_zoom_range);
    matlab2tikz(sprintf("%s/parts/ergebnisse/signal_processing/signal_lp_zoomed.tex", documentation_path), "showInfo", false, 'height', '\fheight', 'width', '\fwidth');
  else
    error("Not implemented yet");
  endif
endfunction

create_water_fit_figures
% create_water_fit("WaterFit", "WaterFit/Index.csv");
% figure(8);
% xlim([28, 40]);
% ylim([3050, 3350]);

% p_an(16);
% figure(1);
% xlim([4.108e-05, 4.125e-05]);
% ylim([0.0132, 0.0142]);
display_single_signal();

% create_peak_combination_figure
% create_signal_figure

% show_surface_plot_probabilities("WaterFit", 21, 20, mean([37.22, 37.32]));
% p_an(3);

% p_an(10);

% for i = 3:28
% printf("Calculating %d\n", i);
%  p_an(i);
%  input(sprintf("Currently displaying %d.csv\n", i));
% endfor
