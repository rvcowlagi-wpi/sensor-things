clear variables; close all; clc;

figure('units', 'normalized', 'OuterPosition', [0.05 0.05 0.8 0.9]);
bar([0 1], [0.5 0.5], 'EdgeColor','white', 'FaceColor','black', ...
	'FaceAlpha',0.6, 'BarWidth',0.1)
xlabel('Outcome $x$', 'FontName', 'Times New Roman', ...
	'FontSize', 20, 'FontWeight', 'bold', 'interpreter', 'latex');
ylabel('$F_X(x)$', 'FontName', 'Times New Roman', ...
	'FontSize', 20, 'FontWeight', 'bold', 'interpreter', 'latex');
ylim([0 1])
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 20;

exportgraphics(ax, 'coin-toss.png', 'Resolution', 300);


figure('units', 'normalized', 'OuterPosition', [0.05 0.05 0.8 0.9]);
bar(1:6, (1/6)*ones(1,6), 'EdgeColor','white', 'FaceColor','black', ...
	'FaceAlpha', 0.6, 'BarWidth',0.1)
xlabel('Outcome $x$', 'FontName', 'Times New Roman', ...
	'FontSize', 20, 'FontWeight', 'bold', 'interpreter', 'latex');
ylabel('$F_X(x)$', 'FontName', 'Times New Roman', ...
	'FontSize', 20, 'FontWeight', 'bold', 'interpreter', 'latex');
ylim([0 1])
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 20;

exportgraphics(ax, 'die-toss.png', 'Resolution', 300);
