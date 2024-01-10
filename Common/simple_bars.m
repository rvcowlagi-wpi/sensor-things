%% 
clear variables; close all; clc;

fontSize_ = 50;
figure('units', 'normalized', 'OuterPosition', [0.05 0.05 0.8 0.9]);
bar([0 1], [0.5 0.5], 'EdgeColor','white', 'FaceColor','black', ...
	'FaceAlpha',0.6, 'BarWidth',0.1)
xlabel('Outcome $x$', 'FontName', 'Times New Roman', ...
	'FontSize', fontSize_, 'FontWeight', 'bold', 'interpreter', 'latex');
ylabel('$F_X(x)$', 'FontName', 'Times New Roman', ...
	'FontSize', fontSize_, 'FontWeight', 'bold', 'interpreter', 'latex');
ylim([0 1])
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = fontSize_;

exportgraphics(ax, 'coin-toss.png', 'Resolution', 300);

%%
clear variables; close all; clc;

fontSize_ = 50;

figure('units', 'normalized', 'OuterPosition', [0.05 0.05 0.8 0.9]);
bar(1:6, (1/6)*ones(1,6), 'EdgeColor','white', 'FaceColor','black', ...
	'FaceAlpha', 0.6, 'BarWidth',0.3)
xlabel('Outcome $x$', 'FontName', 'Times New Roman', ...
	'FontSize', fontSize_, 'FontWeight', 'bold', 'interpreter', 'latex');
ylabel('$F_X(x)$', 'FontName', 'Times New Roman', ...
	'FontSize', fontSize_, 'FontWeight', 'bold', 'interpreter', 'latex');
ylim([0 1])
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = fontSize_;

exportgraphics(ax, 'die-toss.png', 'Resolution', 300);

%%
clear variables; close all; clc;

fontSize_ = 50;
figure('units', 'normalized', 'OuterPosition', [0.05 0.05 0.8 0.9]);
bar(350:10:410, [0.05 0.1 0.2 0.3 0.2 0.1 0.05], ...
	'EdgeColor','white', 'FaceColor','black', ...
	'FaceAlpha', 0.6, 'BarWidth', 0.3)
xlabel('Outcome $x$', 'FontName', 'Times New Roman', ...
	'FontSize', fontSize_, 'FontWeight', 'bold', 'interpreter', 'latex');
ylabel('$F_X(x)$', 'FontName', 'Times New Roman', ...
	'FontSize', fontSize_, 'FontWeight', 'bold', 'interpreter', 'latex');
ylim([0 0.4])
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = fontSize_;

exportgraphics(ax, 'pmf-bean1.png', 'Resolution', 300);


%%
clear variables; close all; clc;

fontSize_ = 50;
figure('units', 'normalized', 'OuterPosition', [0.05 0.05 0.8 0.9]);
bar(350:10:410, [0.1 0.15 0.15 0.2 0.15 0.15 0.1], ...
	'EdgeColor','white', 'FaceColor','black', ...
	'FaceAlpha', 0.6, 'BarWidth', 0.3)
xlabel('Outcome $x$', 'FontName', 'Times New Roman', ...
	'FontSize', fontSize_, 'FontWeight', 'bold', 'interpreter', 'latex');
ylabel('$F_X(x)$', 'FontName', 'Times New Roman', ...
	'FontSize', fontSize_, 'FontWeight', 'bold', 'interpreter', 'latex');
ylim([0 0.4])
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = fontSize_;

exportgraphics(ax, 'pmf-bean2.png', 'Resolution', 300);