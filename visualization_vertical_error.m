graphics_toolkit("fltk");

printf("Starting\n");

global noise_uncertainty_floor;
noise_uncertainty_floor = 0.0004;

function [time, voltage] = load_csv_file(path)
  read_file = csvread(path);

  start_time = read_file(2,3);
  increment_time = read_file(2,4);

  time = start_time + increment_time * read_file(3:2:end, 1);
  voltage = read_file(3:2:end, 2);
endfunction

function result = selector(number)
  result = number > 0;
endfunction

function supersampled = super_sample(array, n)
  supersampled = interp1(1:n:length(array)*n, array, 1:length(array)*n);
endfunction

function mlab_graph(name)
  matlab2tikz(sprintf("/home/budschie/Dokumente/HectorCloud/HectorData/SharedData/VerticalErrorVisualization/%s.tex", name), "showInfo", false);
endfunction

function show_area_of_interest_visualization(time, voltage)
  figure(1);
  plot(time, voltage, "black");
  title("Gemessenes Signal");
  xlabel("U in V");
  ylabel("t in s");

  mlab_graph("Signal");

  figure(2);
  global noise_uncertainty_floor;

  plot(time, voltage, "black");
  hold on;
  plot(time, voltage - noise_uncertainty_floor, "red");
  plot(time, voltage + noise_uncertainty_floor, "green");

  legend("Gemessenes Signal", "Minimaler tatsächlicher Wert", "Maximaler tatsächlicher Wert");
  title("Vertikaler Fehler");
  xlabel("U in V");
  ylabel("t in s");

  mlab_graph("ErroredSignal");

  figure(3);
  global noise_uncertainty_floor;

  min_peak = max(voltage - noise_uncertainty_floor);
  plot(time, repmat(min_peak, 1, length(time)), "blue");
  hold on;
  plot(time, voltage, "black");
  plot(time, voltage - noise_uncertainty_floor, "red");
  plot(time, voltage + noise_uncertainty_floor, "green");

  legend("Kleinste Peak-Höhe", "Gemessenes Signal", "Minimaler tatsächlicher Wert", "Maximaler tatsächlicher Wert");
  title("Kleinstmöglicher Wert des größten Signalwertes");
  xlabel("U in V");
  ylabel("t in s");

  mlab_graph("ErroredPeakThresholdSignal");


  figure(4);
  global noise_uncertainty_floor;

  supsamp_time = super_sample(time, 10);
  supsamp_voltage = super_sample(voltage, 10);
  sorted_zero_voltage = find((supsamp_voltage - min_peak + noise_uncertainty_floor + 0.000001) > 0);

  area(supsamp_time(sorted_zero_voltage), supsamp_voltage(sorted_zero_voltage) + noise_uncertainty_floor, min_peak, "FaceColor", [0.9 0.9 0.5]);
  hold on;

  min_peak = max(voltage - noise_uncertainty_floor);
  plot(time, repmat(min_peak, 1, length(time)), "blue");
  plot(time, voltage, "black");
  plot(time, voltage - noise_uncertainty_floor, "red");
  plot(time, voltage + noise_uncertainty_floor, "green");

  legend("Peak-Bereich", "Minimaler tatsächlicher Wert", "Gemessenes Signal", "Maximaler tatsächlicher Wert", "Kleinste Peak-Höhe");
  title("Peak-Bereich");
  xlabel("U in V");
  ylabel("t in s");

  mlab_graph("ErroredPeakAreaSignal");


endfunction

[time, voltage] = load_csv_file("WaterFit/0.csv");
show_area_of_interest_visualization(time, voltage);
