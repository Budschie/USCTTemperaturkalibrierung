
graphics_toolkit("qt");
1;

function [sampled_x, sampled_y] = get_samples(max_frequency, sampling_time_step, sampling_duration, phase)
  last = 1.0 / max_frequency;
  sampled_x = 0:sampling_time_step:sampling_duration;
  sampled_y = sin((min(sampled_x, last) + phase) * 2 * pi * max_frequency) + sin((min(sampled_x, last) + phase) * 2 * pi * max_frequency * 0.9);
end

% Wrong, interpft has to be used
function target_y = whittaker_shannon_interpolation(target_x, samples_y, T)
  element_index = 1;

  target_y = zeros(1, length(target_x));

  for x = target_x
    value = 0;
    for i = 1:length(samples_y)
      value += sinc((x - i * T) / T) * samples_y(i);
    endfor

    target_y(element_index) = value;
    element_index = element_index + 1;
  endfor
end

%% HUGE COMMENT
%% ACTUALLY NEWVERMIND DONT USE PADDING, IT DECREASES ACCURACY; MAKE SURE THAT SAMPLING DURATION IS 1 / max_frequency
sampling_duration = 2.5e-7;
sampling_time_step = 0.04e-6;
phase = 1.0/3.8e6 * 0.15;

interpft_between_steps = 1000;

% TODO: Try padding at tart
[sampled_x, sampled_y] = get_samples(3.8e6, sampling_time_step, sampling_duration, phase);

[value_max_raw, index_max_raw] = max(sampled_y);
printf("Raw max at %.6e\n", sampled_x(index_max_raw));

interpft_x = 0:(sampling_time_step / interpft_between_steps):sampling_duration;
interpft_y = real(interpft(sampled_y, size(interpft_x)(2)));
% interpft_y = whittaker_shannon_interpolation(interpft_x, sampled_y, sampling_duration);

[value_max_interpolated, index_max_interpolated] = max(interpft_y);

printf("Interpolated max at %.6e\n", interpft_x(index_max_interpolated));

figure();
scatter(sampled_x, sampled_y, 'x');
hold on;
scatter(interpft_x, interpft_y, 'o');
hold on;
[ground_truth_x, ground_truth_y] = get_samples(3.8e6, sampling_time_step / interpft_between_steps, sampling_duration, phase);
scatter(ground_truth_x, ground_truth_y, 's');

[value_max_gt, index_max_gt] = max(ground_truth_y);

printf("Ground truth max at %.6e\n\n", ground_truth_x(index_max_gt));
printf("Thats a difference of %.6e of gt to interpft\n", ground_truth_x(index_max_gt) - interpft_x(index_max_interpolated));
printf("Thats a difference of %.6e of gt to nearest sample\n", ground_truth_x(index_max_gt) - sampled_x(index_max_raw));
