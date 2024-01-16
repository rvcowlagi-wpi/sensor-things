%{
Plots Gaussian PDFs
%}

clear variables; close all; clc


x_	= -5:0.01:5;
mu_ = 0;
s1_	= 0.5;
s2_	= 1;
s3_ = 2;

fx1_ = (1 / sqrt(2*pi*s1_^2))*exp( -(x_ - mu_).^2 ./ (2*s1_^2) );
fx2_ = (1 / sqrt(2*pi*s2_^2))*exp( -(x_ - mu_).^2 ./ (2*s2_^2) );
fx3_ = (1 / sqrt(2*pi*s3_^2))*exp( -(x_ - mu_).^2 ./ (2*s3_^2) );


fontSize_ = 40;
figure('units', 'normalized', 'OuterPosition', [0.05 0.05 0.8 0.9]);
plot(x_, fx1_, 'LineWidth', 4, 'Color', 'k'); hold on; grid on; ax = gca;
% ax.ColorOrderIndex = 1;
plot(x_, fx2_, 'LineWidth', 4, 'LineStyle','--', 'Color', 'k')
% ax.ColorOrderIndex = 1;
plot(x_, fx3_, 'LineWidth', 4, 'LineStyle',':', 'Color', 'k')

legend('\quad$\sigma = 0.5$', '\quad$\sigma = 1$', '\quad$\sigma = 2$', 'Interpreter', 'latex')

xlabel('$x$', 'FontName', 'Times New Roman', ...
	'FontSize', fontSize_, 'FontWeight', 'bold', 'interpreter', 'latex');
ylabel('$f_X(x)$', 'FontName', 'Times New Roman', ...
	'FontSize', fontSize_, 'FontWeight', 'bold', 'interpreter', 'latex');
ylim([0 1])
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = fontSize_;

exportgraphics(ax, 'gaussian-pdf.png', 'Resolution', 300);
