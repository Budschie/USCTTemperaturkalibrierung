% usct3.0_geometry_v1.83_KITnumbering_TAS_rotation_rand_optimized.mat
graphics_toolkit("qt");

global export_mode;
global export_line_width;
% doc, presentation
export_mode="presentation";

global export_path;

if strcmp(export_mode, "doc")
  export_path="/home/budschie/Dokumente/HectorCloud/HectorData/Abschlussdokumentation/parts";
  export_line_width=3;
else
  export_path="/home/budschie/Dokumente/HectorCloud/HectorData/FinalerVortragII/Images/TemperatureDifference";
endif

function save_fig_named(name)
  global export_mode;
  global export_path;

  if strcmp(export_mode, "doc")
    matlab2tikz(sprintf("%s/temperature_analysis/%s.tex", export_path, name), "showInfo", false, 'height', '\fheight', 'width', '\fwidth');
  elseif strcmp(export_mode, "presentation")
    savefig(sprintf("%s/%s", export_path, name));
    print(sprintf("%s/%s.svg", export_path, name), "-dsvg", "-S1920,1080");
  endif
endfunction

set(0, "defaultlinelinewidth", export_line_width);

pkg load statistics;

load("/home/budschie/Dokumente/NextcloudTMG/UltraschallCT/USCTAnalyze/3D-USCT-III-access-script/usct3.0_geometry_v1.83_KITnumbering_TAS_rotation_rand_optimized.mat");
load("/home/budschie/Development/Octave/Hector/Data/_exp878 20.7percent alcohol aver128 (average weak)/TASTemp.mat");

insane_temperatures_x = TASTemperature(2, :);
insane_temperatures_y = TASTemperature(1, :);

[insane_temperatures_sorted, insane_temperatures_sorted_index] = sort(insane_temperatures_y);
removal_radius = 2;

sane_temperatures_sorted_index = insane_temperatures_sorted_index(removal_radius:end-removal_radius);
sane_temperatures_sorted = insane_temperatures_sorted(removal_radius:end-removal_radius);

figure(1);
clf;
plot(1:length(sane_temperatures_sorted), sane_temperatures_sorted);
title("Temperatures measured at TAS");

figure(2);

raw_positions = [TASElements.TASTemperatureSensor];
positions = zeros(length(sane_temperatures_sorted), 3);

for i = 1:length(sane_temperatures_sorted)
  positions(i,:) = raw_positions(:,sane_temperatures_sorted_index(i));
endfor

figure(2);

if strcmp(export_mode, "doc")
  export_radius = 10;
elseif strcmp(export_mode, "presentation")
  export_radius = 160;
endif

scatter3(positions(:,1), positions(:,2), positions(:,3), export_radius, sane_temperatures_sorted, 'filled');
cbar=colorbar();
set(gca, "linewidth", export_line_width, "fontsize", 16);

if strcmp(export_mode, "doc")
  title("XYZ plot of temperature in \\textdegree C against position");
elseif strcmp(export_mode, "presentation")
  xlabel("Breite in m");
  ylabel("Tiefe in m");
  zlabel("Höhe in m");
  title("XYZ-Plot der Temperatur in\\degC");
  title(cbar, "T in \\degC", "fontsize", 22);
endif

save_fig_named("3DTemperature");

figure(3);
scatter(positions(:,3), sane_temperatures_sorted, 'filled');
title("Temperature in °C against Z (up) coordinate");
set(gca, "linewidth", export_line_width, "fontsize", 22);

save_fig_named("ScatterTemperature");

figure(4);
a = boxplot(sane_temperatures_sorted, "linewidth", 10);
title("Boxplot of temperature data");
xlim([0.4, 1.6]);

if strcmp(export_mode, "doc")
  ylabel("$T$ in \\textdegree C");
elseif strcmp(export_mode, "presentation")
  ylabel("T in °C");
endif

set(gca, "linewidth", export_line_width, "fontsize", 22);

save_fig_named("BoxplotTemperature");

printf("Max temperature difference of %.3f°C\n", abs(sane_temperatures_sorted(1) - sane_temperatures_sorted(end)));
