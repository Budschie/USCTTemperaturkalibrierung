% YES
calculate_functions

function use_case_object = new_use_case_short(min_temp, max_temp, name, optimal_wp)
  use_case_object = new_use_case(min_temp, max_temp, name, optimal_wp, 5, 1, 0.0003, 0.999620765, 1.0 - 0.999620765, 20);
endfunction

function use_case_object = new_use_case(min_temp, max_temp, name, optimal_wp, n_water, n_ethanol, linearity, xi_e, k, absolute_mass)
  use_case_object.min_temp = min_temp;
  use_case_object.max_temp = max_temp;
  use_case_object.name = name;
  use_case_object.optimal_wp = optimal_wp;

  % Everything below in SI
  use_case_object.n_water = n_water;
  use_case_object.n_ethanol = n_ethanol;
  use_case_object.linearity = linearity;
  use_case_object.xi_e = xi_e;
  use_case_object.k = k;
  use_case_object.absolute_mass = absolute_mass;
endfunction

global cache_path;
global use_cases;
global export_mode;
global export_line_width;
use_cases = [ new_use_case_short(19, 21, "MK", 15.744), new_use_case_short(18, 28, "ME", 15.287), new_use_case_short(35.2, 37.2, "MR", 12.616), new_use_case_short(18, 45, "MV", 13.747) ];
% doc, presentation
export_mode="presentation";
cache_path=".cache";

global export_path;

if strcmp(export_mode, "doc")
  export_path="/home/budschie/Dokumente/HectorCloud/HectorData/Abschlussdokumentation/parts";
  export_line_width=3;
elseif strcmp(export_mode, "presentation")
  export_path="/home/budschie/Dokumente/HectorCloud/HectorData/FinalerVortragII/Images";
  export_line_width=8;
else
  error("Wat");
endif

% figure(1);
% temperatures = 10:0.05:40;
% wp_1 = 15;
% wp_2 = 12;

% plot(temperatures, calculate_speed_of_sound(temperatures, repmat(wp_1, length(temperatures), 1)), "linewidth", 4);
% hold on;
% plot(temperatures, calculate_speed_of_sound(temperatures, repmat(wp_2, length(temperatures), 1)), "linewidth", 4);
% hold off;

% set(gca, "linewidth", 6, "fontsize", 22);

% title("c gegen T für verschiedene 𝛏");
% legend(sprintf("𝛏=%.2f%%", wp_1), sprintf("𝛏=%.2f%%", wp_2));
% xlabel("T in °C");
% ylabel("c in m/s");

function plot_vertical_lines(x, y, height)
  real_x = repmat(x, 2, 1);
  real_y = zeros(2, length(y));
  real_y(1,:) = y;
  real_y(2,:) = y + height;
  plot(real_x, real_y, "color", "green");
endfunction

function standard_error = get_standard_error(expected, estimate)
  mse_single = (expected - estimate) .* (expected - estimate);
  mse_averaged = mean(mse_single, 2);

  standard_error = sqrt(mse_averaged);
endfunction

% Assumption: Only two crossings
function [first_crossed, second_crossed] = find_crossings(dataset, crossed_number)
  % Calculate whether we are above or below our crossed number
  below_number = dataset < crossed_number;
  % Find the largest and smallest index, for which below_number == 1

  first_crossed = -1;
  second_crossed = -1;

  for i=1:length(below_number)
    if below_number(i) == 1
      if first_crossed == -1
        first_crossed = i;
      endif

      second_crossed = i;
    endif
  endfor
endfunction

% Function maps wp to minimum required mixing accuracy on a conceptual level
function [wp, minimum_required_mixing_accuracies] = calculate_min_required_mixing_accuracies(temperatures, wp_min, wp_max, wp_deviation_min, wp_deviation_max, wp_deviation_resolution, wp_error_skips, crossing)
  % [temperature_meshgrid, wp_meshgrid] = meshgrid(temperature_range, wp_range);

  % 1nd index, 2nd : is wp-coherent

  % baselines = calculate_bulk_speed_of_sound(temperatures, wp_meshgrid);

  % minimum_required_mixing_accuracies = zeros(length(wp_range), 1);

  all_required_wp = wp_min-wp_deviation_min:wp_deviation_resolution:wp_max+wp_deviation_max;

  [temp_grid, wp_grid] = meshgrid(temperatures, all_required_wp);
  all_required_sos = calculate_bulk_speed_of_sound(temp_grid, wp_grid);

  min_max_buffer_size = wp_deviation_min/wp_deviation_resolution;
  start_real_index = min_max_buffer_size + 1;
  stop_real_index = length(all_required_wp) - min_max_buffer_size;

  baseline_indices=start_real_index:wp_error_skips:stop_real_index;

  wp = all_required_wp(baseline_indices);
  minimum_required_mixing_accuracies = zeros(length(wp), 1);

  deviation_coordinate_space = -wp_deviation_min:wp_deviation_resolution:wp_deviation_max;

  for j=1:length(baseline_indices)
    i = baseline_indices(j);
    current_wp = wp(j);
    baseline_repeated = repmat(all_required_sos(i,:), 2 * min_max_buffer_size + 1, 1);

    errors = get_standard_error(baseline_repeated, all_required_sos(i-min_max_buffer_size:i+min_max_buffer_size,:));
    [first_crossed, second_crossed] = find_crossings(errors, crossing);
    crossed_values = [deviation_coordinate_space(first_crossed), all_required_wp(second_crossed)];
    abs_min_crossed_value = min(abs(crossed_values));

    minimum_required_mixing_accuracies(j) = abs_min_crossed_value;
    % deviated_wps = current_wp + wp_deviation;

    % [temperature_2_meshgrid, wp_2_meshgrid] = meshgrid(temperature_range, deviated_wps);
    % estimates = calculate_bulk_speed_of_sound(temperature_2_meshgrid, wp_2_meshgrid);

    % current_baseline_repeated = repmat(baselines(i,:), length(deviated_wps), 1);

    % errors = get_standard_error(current_baseline_repeated, estimates);

    % [first_crossed, second_crossed] = find_crossings(errors, 1);
    % crossed_values = [wp_deviation(first_crossed), wp_deviation(second_crossed)];
    % abs_min_crossed_value = min(abs(crossed_values));

    % minimum_required_mixing_accuracies(i) = abs_min_crossed_value;
  endfor
endfunction

function save_visualize_mean_abs_diff(name)
  print(sprintf("/home/budschie/Dokumente/HectorCloud/HectorData/FinalerVortragII/Images/TempDepNew/%s.svg", name), "-S1920,1080");
endfunction

function visualize_mean_absolute_difference(data_x, data_y, mu, error_skips)
  clf;
  set(gca, "fontsize", 18);
  hold on;
  title("c gegen T für Beispiel-Ethanol-Wasser-Mischung", "fontsize", 40);
  xlabel("T in \\degC");
  ylabel("c in m/s");
  plot(data_x, data_y, "linewidth", 4);
  save_visualize_mean_abs_diff("1full_plot");
  clf;
  set(gca, "fontsize", 18);
  hold on;
  title("c gegen T für Beispiel-Ethanol-Wasser-Mischung");
  xlabel("T in \\degC");
  ylabel("c in m/s");

  skip_indices=1:error_skips:length(data_x);
  diff = data_y - mu;
  % helper_line = mu + diff / 2;

  plot(data_x(skip_indices), data_y(skip_indices), "linewidth", 4, "marker", "o", "linestyle", "none");
  save_visualize_mean_abs_diff("2reduced_plot");

  % plot(data_x, data_y, "linewidth", 4, "linestyle", "--");
  plot(data_x, repmat(mu, length(data_x), 1), "linewidth", 4);
  save_visualize_mean_abs_diff("3median");

  plot_vertical_lines(data_x(1:error_skips:length(data_x)), data_y(1:error_skips:length(data_x)), -diff(1:error_skips:length(data_x)));  % hold on;
  save_visualize_mean_abs_diff("4diff");
  % test = errorbar(data_x, helper_line, repmat(0, length(diff), 1), diff / 2, "#~>");

  hold off;
endfunction

function visualize_difference(data_x, data_y, first_index, second_index, title_text, text_shift_factor, figure_indices)
  figure(figure_indices(1));
  clf;
  hold on;
  plot(data_x, data_y, "linewidth", 4);
  set(gca, "linewidth", 6);
  title(title_text);
  xlabel("T in \\degC");
  ylabel("c in m/s");

  figure(figure_indices(2));
  clf;
  plot(data_x, data_y, "linewidth", 4);
  hold on;
  plot([data_x(first_index), data_x(second_index)], [data_y(first_index), data_y(second_index)], "linewidth", 13, "linestyle", "none", "marker", "o", "color", "red");

  xlabel("T in \\degC");
  ylabel("c in m/s");
  set(gca, "linewidth", 6);
  title(title_text);

  figure(figure_indices(3));
  clf;
  plot(data_x, data_y, "linewidth", 4);
  hold on;
  plot([data_x(first_index), data_x(second_index)], [data_y(first_index), data_y(second_index)], "linewidth", 13, "linestyle", "none", "marker", "o", "color", "red");

  % Left horizontal line
  first_bound = data_y(first_index);
  second_bound = data_y(second_index);

  plot([data_x(1), data_x(end)], [first_bound, first_bound], "linewidth", 3, "color", "#F0D752");
  % Right horizontal line
  plot([data_x(1), data_x(end)], [second_bound, second_bound], "linewidth", 3, "color", "#F0D752");
  % Middle vertical line
  plot([data_x(floor(end/2)), data_x(floor(end/2))], [first_bound, second_bound], "linewidth", 3, "color", "#F0D752", "marker", "o", "linestyle", "--");

  % Text right to that vertical line
  diff_string = sprintf("\\Deltac = %.3f m/s", abs(second_bound - first_bound));
  text(data_x(floor(end/2)) + (data_x(end) - data_x(1)) * text_shift_factor, (first_bound + second_bound) / 2, diff_string, "fontsize", 17, "horizontalalignment", "center");

  xlabel("T in \\degC");
  ylabel("c in m/s");
  set(gca, "linewidth", 6);
  title(title_text);
endfunction

function complete_boundary_difference_visualization
  wp = 14;
  T_1 = 26:0.1:30;
  T_2 = 20:0.1:30;

  sos_1 = calculate_speed_of_sound(wp, T_1);
  sos_2 = calculate_speed_of_sound(wp, T_2);

  visualize_difference(T_1, sos_1, 1, length(T_1), "Differenz der Randwerte", -0.24, [3, 4, 5]);
  visualize_difference(T_2, sos_2, 1, length(T_2), "Differenz der Randwerte", -0.24, [6, 7, 8]);

  path="/home/budschie/Dokumente/HectorCloud/HectorData/FinalerVortragII/Images/BoundaryDifference";

  for i=[3:8]
    figure(i);
    savefig(sprintf("%s/Figure%d", path, i));
    print(sprintf("%s/Figure%d.svg", path, i), "-dsvg", "-S1920,1080");
  endfor
endfunction

function complete_min_max_difference_visualization
  wp = 14;
  T_2 = 20:0.1:30;

  sos_2 = calculate_speed_of_sound(wp, T_2);

  [min_val, min_index] = min(sos_2);
  [max_val, max_index] = max(sos_2);

  visualize_difference(T_2, sos_2, min_index, max_index, "Differenz der Extrema", -0.24, [9, 10, 11]);

  path="/home/budschie/Dokumente/HectorCloud/HectorData/FinalerVortragII/Images/MinMaxDifference";

  for i=[9:11]
    figure(i);
    savefig(sprintf("%s/Figure%d", path, i));
    print(sprintf("%s/Figure%d.svg", path, i), "-dsvg", "-S1920,1080");
  endfor

endfunction

function calculate_min_required_mixing_accuracy_graph(crossing)
  global use_cases;
  global cache_path;

  for i=1:length(use_cases)
      figure(12);
      clf;

      temperature_range = use_cases(i).min_temp:0.1:use_cases(i).max_temp;

      [wp,min_required_mixing_accuracy] = calculate_min_required_mixing_accuracies(temperature_range, 4, 16, 4, 4, 0.00005, round(0.1/0.00005), crossing);
      plot(wp,min_required_mixing_accuracy);

      savefig(sprintf("%s/min_required_accuracy_%.2f_%s.fig", cache_path, crossing, use_cases(i).name));
  endfor
endfunction

function delta_xi = calc_wrong_raw_ethanol_wp(k, xi_e, xi)
  delta_xi = xi .* (k ./ xi_e);
endfunction

function [lin_err, min_eth_wp, max_eth_wp] = calc_linearity_error(mass, xi, n_water, n_ethanol, linearity)
  n_dot_lin_wat = linearity .* n_water;
  n_dot_lin_eth = linearity .* n_ethanol;

  mass_eth = xi .* mass;
  mass_wat = (1 - xi) .* mass;

  min_eth_wp = (mass_eth - n_dot_lin_eth) ./ (mass_eth - n_dot_lin_eth + mass_wat + n_dot_lin_wat);
  max_eth_wp = (mass_eth + n_dot_lin_eth) ./ (mass_eth + n_dot_lin_eth + mass_wat - n_dot_lin_wat);

  diff_min = abs(min_eth_wp - xi);
  diff_max = abs(max_eth_wp - xi);

  lin_err = (diff_min > diff_max) .* diff_min + (diff_min <= diff_max) .* diff_max;
endfunction

% Maybe also add prob distribution
function [full_error, error_without_xi] = calculate_full_error(optimal_wp, real_wp, t_min, t_max)
  t_area = t_min:0.001:t_max;
  median_c = median(calculate_bulk_speed_of_sound(t_area, repmat(optimal_wp, length(t_area), 1)));
  full_error = mean(abs(calculate_bulk_speed_of_sound(t_area, repmat(real_wp, length(t_area), 1)) - median_c));
  error_without_xi = mean(abs(calculate_bulk_speed_of_sound(t_area, repmat(optimal_wp, length(t_area), 1)) - median_c));
endfunction

function new_figure = merge_figures(figure_list, color_list, thickness_list, axes_closure)
  new_figure = figure();

  for i = 1:length(figure_list)
    current = figure_list(i);
    axis = findobj(current, "Type", "Axes");

    set(get(axis, "Children")(1), "Color", strvcat(color_list(i)));
    set(get(axis, "Children")(1), "linewidth", thickness_list(i));
    % set(get(axis, "Children")(1), "linestyle", thickness_list(i));
    axes_closure(axis);

    if i==1
      figure_axis = copyobj(axis, new_figure);
    else
      copyobj(get(axis, "Children"), figure_axis);
    endif

    close(current);
  endfor
endfunction

function std_err = calculate_mixing_error(assumed_wp, real_wp, min_temp, max_temp)
  temps = min_temp:0.02:max_temp;
  std_err = get_standard_error(calculate_bulk_speed_of_sound(temps, repmat(real_wp * 100, length(temps), 1))', calculate_bulk_speed_of_sound(temps, repmat(assumed_wp * 100, length(temps), 1))');
endfunction

function display_mixing_error()
  global use_cases;

  for i=1:length(use_cases)
     ethanol_wp_error_pos = calc_wrong_raw_ethanol_wp(use_cases(i).k, use_cases(i).xi_e, use_cases(i).optimal_wp / 100);
     ethanol_wp_error_neg = calc_wrong_raw_ethanol_wp(-use_cases(i).k, use_cases(i).xi_e, use_cases(i).optimal_wp / 100);

     lin_err_min = 1000;
     lin_err_max = 0;

     errors = [ethanol_wp_error_neg, ethanol_wp_error_pos];

       for j = 1:length(errors)
        [linearity_error_temp, lin_min_tmp, lin_max_tmp] = calc_linearity_error(use_cases(i).absolute_mass,
          use_cases(i).optimal_wp / 100 + errors(j), ...
          use_cases(i).n_water, ...
         use_cases(i).n_ethanol, ...
          use_cases(i).linearity);
          % lin_err = max(linearity_error_temp, lin_err)(1);
          lin_err_min = min(lin_err_min, lin_min_tmp);
          lin_err_max = max(lin_err_max, lin_max_tmp);
       endfor

       [std_err_min, optimal_error] = ...
          calculate_full_error(use_cases(i).optimal_wp, ...
          lin_err_min * 100, ...
          use_cases(i).min_temp, ...
          use_cases(i).max_temp);

       std_err_max = ...
          calculate_full_error(use_cases(i).optimal_wp, ...
          lin_err_max * 100, ...
          use_cases(i).min_temp, ...
          use_cases(i).max_temp);

       printf("use case with %.1f to %.1f: error of %.3e compared to %.3e\n", use_cases(i).min_temp, use_cases(i).max_temp, max(std_err_min, std_err_max), optimal_error);
  endfor
endfunction

function calculate_mixing_errors_values
  global use_cases;

  for i=1:length(use_cases)
    current = use_cases(i);

     ethanol_wp_error_pos = calc_wrong_raw_ethanol_wp(use_cases(i).k, use_cases(i).xi_e, use_cases(i).optimal_wp / 100);
     ethanol_wp_error_neg = calc_wrong_raw_ethanol_wp(-use_cases(i).k, use_cases(i).xi_e, use_cases(i).optimal_wp / 100);

     lin_err_min = 1000;
     lin_err_max = 0;

     errors = [ethanol_wp_error_neg, ethanol_wp_error_pos];

     [linearity_error, min_wp, max_wp] = calc_linearity_error(use_cases(i).absolute_mass,
       use_cases(i).optimal_wp / 100, ...
       use_cases(i).n_water, ...
       use_cases(i).n_ethanol, ...
       use_cases(i).linearity);

       for j = 1:length(errors)
        [linearity_error_temp, lin_min_tmp, lin_max_tmp] = calc_linearity_error(use_cases(i).absolute_mass,
          use_cases(i).optimal_wp / 100 + errors(j), ...
          use_cases(i).n_water, ...
         use_cases(i).n_ethanol, ...
          use_cases(i).linearity);
          % lin_err = max(linearity_error_temp, lin_err)(1);
          lin_err_min = min(lin_err_min, lin_min_tmp);
          lin_err_max = max(lin_err_max, lin_max_tmp);
       endfor

       standalone_wp_error = max(calculate_mixing_error(current.optimal_wp / 100, current.optimal_wp / 100 + errors(1), current.min_temp, current.max_temp),...
       calculate_mixing_error(current.optimal_wp / 100, current.optimal_wp / 100 + errors(2), current.min_temp, current.max_temp));

       standalone_linearity_error = max(calculate_mixing_error(current.optimal_wp / 100, min_wp, current.min_temp, current.max_temp),...
       calculate_mixing_error(current.optimal_wp / 100, max_wp, current.min_temp, current.max_temp));

       combined_error = max(calculate_mixing_error(current.optimal_wp / 100, lin_err_min, current.min_temp, current.max_temp),...
       calculate_mixing_error(current.optimal_wp / 100, lin_err_max, current.min_temp, current.max_temp));


       printf("Use case %s:\n- Standalone  wp err: %.5e\n- Standalone lin err: %.5e\n- Combined       err: %.5e\n\n", current.name, standalone_wp_error, standalone_linearity_error, combined_error);
  endfor
endfunction

function linestyle_change_closure(axes)
  child = get(axes, "Children");
  set(child, "linestyle", "--");

  mid_x = get(child, "xdata")(floor(end/2));
  mid_y = get(child, "ydata")(floor(end/2));

  % text(mid_x, mid_y, "Test", "fontsize", 30);
endfunction

function visualize_min_required_accuracy_graph(crossings)
  global use_cases;
  global cache_path;
  global export_mode;
  global export_path;
  global export_line_width;

  % min_required_accuracy_name="MinRequiredAccuracy";

  for i=1:length(use_cases)
      figure(12);
      clf;

      crossings_indices = zeros(length(crossings), 1);

      for j=1:length(crossings)
        crossings_indices(j) = openfig(sprintf("%s/min_required_accuracy_%.2f_%s.fig", cache_path, crossings(j), use_cases(i).name));
      endfor

      main_plot = merge_figures(crossings_indices, cellstr(["#d95319"; "#edb120"]), repmat(export_line_width*0.5, length(crossings), 1), @linestyle_change_closure);

      wp=4:0.01:16;

      neg_xi_err = calc_wrong_raw_ethanol_wp(...
            repmat(-use_cases(i).k, length(wp), 1),...
            repmat(use_cases(i).xi_e, length(wp), 1),...
            wp' ./ 100);

      pos_xi_err = calc_wrong_raw_ethanol_wp(...
            repmat(use_cases(i).k, length(wp), 1),...
            repmat(use_cases(i).xi_e, length(wp), 1),...
            wp' ./ 100);

       errs = [neg_xi_err, pos_xi_err];

       linearity_error = zeros(length(pos_xi_err), 1);

       for l=1:size(errs)(2)
        [_, lin_min_tmp, lin_max_tmp] = calc_linearity_error( ...
        repmat(use_cases(i).absolute_mass, length(wp), 1), ...
        wp' ./ 100 + errs(:,l), ...
        repmat(use_cases(i).n_water, length(wp), 1), ...
        repmat(use_cases(i).n_ethanol, length(wp), 1), ...
        repmat(use_cases(i).linearity, length(wp), 1));

        linearity_error_to_max = max(abs(lin_min_tmp - wp' ./ 100), abs(lin_max_tmp - wp' ./ 100));

        linearity_error = max(linearity_error, linearity_error_to_max);
       endfor


        % delta_xi_error = max(
        %   calc_wrong_raw_ethanol_wp(
        %     repmat(-use_cases(i).k, length(wp), 1),
        %     repmat(use_cases(i).xi_e, length(wp), 1),
        %   wp' ./ 100));

      hold on;

      % set(gca,'yscale','log');
      % set(gca,'xscale','log');

      legend("$\\sigma = 0.1$ m/s");

      grid on;

      if strcmp(export_mode, "doc")
        % xlim([4, 16]);
        title(sprintf("Minimal required mixing accuracy\\\\ for $%.2f$ \\textdegree C $< T < %.2f$ \\textdegree C", use_cases(i).min_temp, use_cases(i).max_temp));
        xlabel("$\\xi_{assumed}$ in w/w \\%");
        ylabel("$\\Delta \\xi$ in w/w \\%");
        matlab2tikz(sprintf("%s/material_und_methoden/c_errors/min_required_accuracy_%s.tex", export_path, use_cases(i).name), "showInfo", false, 'height', '\fheight', 'width', '\fwidth');
      else
        legend("\\sigma = 0.1 m/s");
        title(sprintf("Minimal benötigte Mischgenauigkeit für %.2f°C < T < %.2f°C", use_cases(i).min_temp, use_cases(i).max_temp));
        xlabel("\\xi_a_s_s_u_m_e_d in w/w %");
        ylabel("\\Delta\\xi in w/w %");
        print(sprintf("%s/MixingAccuracy/MinRequiredAccuracy.svg", export_path), "-dsvg", "-S1920,1080");
      endif


      set(gca, "yscale", "log");

      if strcmp(export_mode, "doc")
        plot(wp, max(abs(neg_xi_err), abs(pos_xi_err)) * 100, "blue");
        plot(wp, linearity_error * 100, "magenta");
        % xlim([4, 16]);
        title(sprintf("Minimal required mixing accuracy\\\\ for $%.2f$ \\textdegree C $< T < %.2f$ \\textdegree C", use_cases(i).min_temp, use_cases(i).max_temp));
        xlabel("$\\xi_{assumed}$ in w/w \\%");
        ylabel("$\\Delta \\xi$ in w/w \\%");
        legend("$\\sigma = 1$ m/s", "$\\sigma = 0.1$ m/s", "$\\xi_e$ errors", "Scale and $\\xi_e$ errors");
        matlab2tikz(sprintf("%s/ergebnisse/c_errors/min_required_accuracy_%s.tex", export_path, use_cases(i).name), "showInfo", false, 'height', '\fheight', 'width', '\fwidth');
      else
        eth_q_plot = plot(wp, max(abs(neg_xi_err), abs(pos_xi_err)) * 100, "blue");

        title(sprintf("Minimal benötigte Mischgenauigkeit für %.2f°C < T < %.2f°C", use_cases(i).min_temp, use_cases(i).max_temp));
        xlabel("\\xi_a_s_s_u_m_e_d in w/w %");
        ylabel("\\Delta\\xi in w/w %");
        legend("\\sigma = 0.1 m/s", "Fehler durch Ethanolqualität");
        print(sprintf("%s/MixingAccuracy/MinRequiredAccuracyErrors.svg", export_path), "-dsvg", "-S1920,1080");
        delete(eth_q_plot);
        plot(wp, linearity_error * 100, "magenta");
        legend("\\sigma = 0.1 m/s", "Fehler durch Ethanolqualität und Waage");
        print(sprintf("%s/MixingAccuracy/MinRequiredAccuracyErrors2.svg", export_path), "-dsvg", "-S1920,1080");
      endif

      % xlabel("$\xi_{assumed}$ in w/w");
      % ylabel("$\Delta \xi$ in w/w");
      grid on;

      set(gca, "linewidth", export_line_width, "fontsize", 22);
      set(0, "defaultlinelinewidth", export_line_width);

      grid off;

      % savefig(sprintf("%s/%s", path, min_required_accuracy_name));
      % print(sprintf("%s/%s.svg", path, min_required_accuracy_name), "-dsvg", "-S1920,1080");
      % matlab2tikz("min_required_accuracy.tex", "showInfo", false, 'height', '\fheight', 'width', '\fwidth');
  endfor
endfunction

function diff = xi_plus_diff(xi, xi_e, k)
  % mr = xi ./ xi_e;
  % m_wr = 1 - mr;

  % m_e = (xi_e + k) .* mr;
  % xi_plus = (xi_e + k) .* mr;

  % xi_plus_2 = xi .* (k ./ xi_e);

  % xi_test = (xi_e * mr) / (m_wr + mr);
  % diff = xi - xi_plus;
  diff = xi .* k ./ xi_e;
endfunction

function test_diffs_3d
  xi = 0.05:0.001:0.2;

  changing_k=-0.05:0.0005:0.05;
  changing_xi=0.05:0.0005:0.95;
  [mg_xi, mg_k] = meshgrid(changing_xi, changing_k);
  sizes_mg = size(mg_xi);

  for i=1:length(xi)
    rep_xi = repmat(xi, sizes_mg(1), sizes_mg(2));

    differences = xi_plus_diff(rep_xi, mg_xi, mg_k);


  endfor
endfunction

function test_diffs_graph_changing_xi
  changing_xi = 0.05:0.001:0.2;

  k=0.05;
  xi_e=0.95;

  rep_k = repmat(xi, length(changing_xi), 1);
  rep_xi_e = repmat(xi_e, length(changing_xi), 1);

  differences = xi_plus_diff(changing_xi', rep_xi_e, rep_k);

  figure(1003);
  clf;
  xlabel("xi_eth");
  ylabel("k");
  plot(changing_xi, differences);
  shading interp;
endfunction


function test_diffs_graph_changing_xi_e
  xi = 0.12;

  k=-0.05;
  changing_xi_e=0.05:0.001:0.95;

  rep_xi = repmat(xi, length(changing_xi_e), 1);
  rep_k = repmat(k, length(changing_xi_e), 1);

  differences = xi_plus_diff(rep_xi, changing_xi_e', rep_k);

  figure(1005);
  clf;
  xlabel("xi_eth");
  ylabel("k");
  plot(changing_xi_e, differences);
  shading interp;
endfunction

function test_diffs_graph_changing_k
  xi = 0.12;

  changing_k=-0.05:0.0005:0.05;
  xi_e=0.95;

  rep_xi = repmat(xi, length(changing_k), 1);
  rep_xi_e = repmat(xi_e, length(changing_k), 1);

  differences = xi_plus_diff(rep_xi, rep_xi_e, changing_k');

  figure(1002);
  clf;
  xlabel("xi_eth");
  ylabel("k");
  plot(changing_k, differences);
  shading interp;
endfunction

function test_diffs
  xi = 0.12;

  changing_k=-0.05:0.0005:0.05;
  changing_xi=0.05:0.0005:0.95;

  [mg_xi, mg_k] = meshgrid(changing_xi, changing_k);
  sizes_mg = size(mg_xi);
  rep_xi = repmat(xi, sizes_mg(1), sizes_mg(2));

  differences = xi_plus_diff(rep_xi, mg_xi, mg_k);

  figure(1001);
  clf;
  xlabel("xi_eth");
  ylabel("k");
  surf(mg_xi, mg_k, abs(differences));
  shading interp;
endfunction

function visualize_single_required_accuracy
  global use_cases;
  global export_mode;
  global export_path;
  global export_line_width;

  % path="/home/budschie/Dokumente/HectorCloud/HectorData/FinalerVortragII/Images/MixingAccuracy";

  for i=1:length(use_cases)
    temperature_range = use_cases(i).min_temp:0.1:use_cases(i).max_temp;

    figure(100);
    clf;
    wp=use_cases(i).optimal_wp;
    wp_deviations = -4:0.001:4;

    [mesh_temp, mesh_wp] = meshgrid(temperature_range, wp_deviations + wp);

    wp_repeated = repmat(wp, length(temperature_range), 1);
    baseline = calculate_bulk_speed_of_sound(temperature_range, wp_repeated);
    sos_deviations = calculate_bulk_speed_of_sound(mesh_temp, mesh_wp);

    errors = get_standard_error(repmat(baseline, 1, length(wp_deviations))', sos_deviations);

    if strcmp(export_mode, "presentation")
      errors_tmp = errors;
      % Slightly offset errors for log plot
      errors = max(errors_tmp, repmat(1e-3, length(errors_tmp), 1));
    endif

    plot(wp_deviations, errors, "linewidth", export_line_width);
    hold on;
    first_x = wp_deviations(1);
    last_x = wp_deviations(end);

    % plot([first_x, last_x], [1, 1], "linewidth", export_line_width/2, "--");
    plot([first_x, last_x], [0.1, 0.1], "linewidth", export_line_width/2, "--");

    if strcmp(export_mode, "presentation")
      title(sprintf("\\sigma_c durch Mischfehler bei \\xi=%.1f %% für %.2f \\degC \\leq T \\leq %.2f \\degC", use_cases(i).optimal_wp, use_cases(i).min_temp, use_cases(i).max_temp));
      xlabel("\\Delta\\xi in w/w");
      ylabel("Standardfehler \\sigma_c in m/s");
      legend("c Abweichungen", "\\sigma_c = 0.1 m/s");

      set(gca, "yscale", "log");
      ylim([1e-2, 5]);
      xlim([-0.5, 0.5]);
      xticks(-3:0.05:3);
      print(sprintf("%s/MixingAccuracy/SingleRequiredAccuracy.svg", export_path), "-dsvg", "-S1920,1080");
    elseif strcmp(export_mode, "doc")
      title(sprintf("$\\sigma_c$ due to mixing error at \\xi=%.1f \\%% \\\\ for $%.2f$ \\textdegree C $\\leq T \\leq %.2f$ \\textdegree C", use_cases(i).optimal_wp, use_cases(i).min_temp, use_cases(i).max_temp));
      xlabel("$\\Delta \\xi$ in w/w");
      ylabel("$\\sigma_c$ in m/s");
      legend("$c$ standard error", "standard error of 0.1 m/s");
      grid off;
      matlab2tikz(sprintf("%s/material_und_methoden/c_errors/wp_error_%s.tex", export_path, use_cases(i).name), "showInfo", false, 'height', '\fheight', 'width', '\fwidth');

      xlim([-0.4,0.4]);
      ylim([0,1.1]);
      grid on;
      matlab2tikz(sprintf("%s/material_und_methoden/c_errors/wp_error_%s_zoomed.tex", export_path, use_cases(i).name), "showInfo", false, 'height', '\fheight', 'width', '\fwidth');
      grid off;
    endif
  endfor

  % savefig(sprintf("%s/SingleRequiredAccuracy", path));
  % print(sprintf("%s/SingleRequiredAccuracy.svg", path), "-dsvg", "-S1920,1080");
endfunction

function visualize_complete_tong_fit
  temps = 20:0.4:40;
  wps = 0:0.2:20;

  [temp_mg, wp_mg] = meshgrid(temps, wps);
  sos = calculate_bulk_speed_of_sound(temp_mg, wp_mg);
  % graphics_toolkit("fltk");
  figure(20);
  clf;
  surf(temp_mg, wp_mg, sos);
  set(gca, "fontsize", 18);
  title("Speed of sound according to J. Tong et al.");
  xlabel("T in \\degC");
  ylabel("\\xi in w/w");
  zlabel("c in m/s");
  shading interp;
  cbar = colorbar;
  ylabel(cbar, "c in m/s");
  set(cbar, "fontsize", 18);

  % savefig("/home/budschie/Dokumente/HectorCloud/HectorData/FinalerVortragII/Images/Literaturrecherche/TongFit");
  print("/home/budschie/Dokumente/HectorCloud/HectorData/FinalerVortragII/Images/Literaturrecherche/TongFit.svg", "-dsvg", "-S1920,1080");
  print("/home/budschie/Dokumente/HectorCloud/HectorData/FinalerVortragII/Images/Literaturrecherche/TongFit.png", "-dpng", "-S1920,1080");
  % matlab2tikz("tong_fit.tex", "showInfo", false, 'height', '\fheight', 'width', '\fwidth');
endfunction

function __visualize_vertical_error_set_axis_and_print(name)
  title("Vertical error");
  xlabel("t in s");
  ylabel("U in V");

  axis([4.5017e-5, 4.5182e-5, 0.0268, 0.0345]);

  print(sprintf("%s/%s.png", "/home/budschie/Dokumente/HectorCloud/HectorData/FinalerVortragII/Images/VerticalError/Theory", name), "-dpng", "-S1920,1080");
endfunction

function visualize_vertical_error(t_1, U_1)
  noise_uncertainty_floor = 0.0004;

  figure(50);
  clf;

  plot(t_1, U_1);
  __visualize_vertical_error_set_axis_and_print("signal");

  figure(51);
  clf;
  hold on;
  plot(t_1, U_1);
  plot(t_1, U_1 - noise_uncertainty_floor, "green");
  __visualize_vertical_error_set_axis_and_print("lower");

  figure(52);
  clf;
  hold on;
  plot(t_1, U_1);
  plot(t_1, U_1 - noise_uncertainty_floor, "green");
  plot(t_1, U_1 + noise_uncertainty_floor, "red");
  __visualize_vertical_error_set_axis_and_print("lower_higher");

  figure(53);
  clf;
  hold on;
  plot(t_1, U_1 - noise_uncertainty_floor, "green");
  plot(t_1, U_1 + noise_uncertainty_floor, "red");
  __visualize_vertical_error_set_axis_and_print("lower_higher_no_signal");

  figure(54);
  clf;
  hold on;
  plot(t_1, U_1 - noise_uncertainty_floor, "green");
  plot(t_1, U_1 + noise_uncertainty_floor, "red");

  lowest_possible_peak_height = max(U_1) - noise_uncertainty_floor;

  plot([t_1(1), t_1(end)], repmat(lowest_possible_peak_height, 2, 1));
    __visualize_vertical_error_set_axis_and_print("lower_higher_area_of_peak");



  hold off;
endfunction

[t_1, U_1] = load_csv_file("WaterFit/21.csv");

% calculate_mixing_errors_values

% test_diffs;

% visualize_min_required_accuracy_graph([1, 0.1]);

% [interpft_t_1, interpft_U_1] = perform_interpft(t_1, U_1, 30);
% visualize_vertical_error(interpft_t_1, interpft_U_1);

% visualize_complete_tong_fit

% visualize_single_required_accuracy;
% min_required_accuracy_graph;
% complete_boundary_difference_visualization;
% complete_min_max_difference_visualization;

% figure(2);

% sos_wp_1 = calculate_speed_of_sound(temperatures, repmat(wp_1, length(temperatures), 1));

% temperatures = 10:0.25:12;
% sos_example_1 = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5];
% sos_example_2 = [1, 1, 1, 1, 1, 1, 1, 5, 5];

% temperature_range=19:0.01:21;
% sos_ex = calculate_bulk_speed_of_sound(temperature_range, repmat(15.744, length(temperature_range), 1));
% visualize_mean_absolute_difference(temperature_range, sos_ex, median(sos_ex), 10);

% set(gca, "linewidth", 6, "fontsize", 22);
% xlabel("T");
% ylabel("c");
% title("Veranschaulichung der Temperaturabhängigkeitsmetriken (erfundene Daten)");
