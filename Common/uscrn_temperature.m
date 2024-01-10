%{
Plots temperature data from USCRN Durham NH station over five years.
%}

clear variables; close all; clc;

data2019	= readtable("CRND0103-2019-NH_Durham_2_N.txt");
data2020	= readtable("CRND0103-2020-NH_Durham_2_N.txt");
data2021	= readtable("CRND0103-2021-NH_Durham_2_N.txt");
data2022	= readtable("CRND0103-2022-NH_Durham_2_N.txt");
data2023	= readtable("CRND0103-2023-NH_Durham_2_N.txt");

temp2019	= table2array(data2019(:,6));
temp2020	= table2array(data2020(:,6));
temp2021	= table2array(data2021(:,6));
temp2022	= table2array(data2022(:,6));
temp2023	= table2array(data2023(:,6));

temp5yr		= [temp2019; temp2020; temp2021; temp2022; temp2023];

temp5yr(temp5yr == -9999) = [];

figure('units', 'normalized', 'OuterPosition', [0.05 0.05 0.8 0.9]);
histogram(temp5yr, 40, 'EdgeColor', 'white', 'FaceColor','black')
xlabel('Temperature ($^\circ$C)', 'FontName', 'Times New Roman', ...
	'FontSize', 20, 'FontWeight', 'bold', 'interpreter', 'latex');
ylabel('No. of occurences', 'FontName', 'Times New Roman', ...
	'FontSize', 20, 'FontWeight', 'bold', 'interpreter', 'latex');
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 20;

exportgraphics(ax, 'uscrn-temp-histogram.png', 'Resolution', 300);

close all;

figure('units', 'normalized', 'OuterPosition', [0.05 0.05 0.8 0.9]);
histogram(temp5yr-273.16, 40, 'EdgeColor', 'black', ...
	'DisplayStyle','stairs', 'Normalization','probability', 'LineWidth', 2)
xlabel('Temperature ($^\circ$K)', 'FontName', 'Times New Roman', ...
	'FontSize', 20, 'FontWeight', 'bold', 'interpreter', 'latex');
ylabel('$f_X(x)$', 'FontName', 'Times New Roman', ...
	'FontSize', 20, 'FontWeight', 'bold', 'interpreter', 'latex');
ylim([0 1])
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 20;

exportgraphics(ax, 'uscrn-temp-pdf.png', 'Resolution', 300);
