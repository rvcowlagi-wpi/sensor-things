%%
clear variables; close all; clc;

fontSize_ = 40;

n_	= 10;
p_	= 0.2;
for k_ = 1:n_
	FXBin(k_)	= nchoosek(n_, k_)* p_^k_ * (1 - p_)^(n_ - k_);
end

figure('units', 'normalized', 'OuterPosition', [0.05 0.05 0.8 0.9]);
bar(1:n_, FXBin, 'EdgeColor','white', 'FaceColor','black', ...
	'FaceAlpha', 0.6, 'BarWidth', 0.3)
xlabel('Outcome $x$', 'FontName', 'Times New Roman', ...
	'FontSize', fontSize_, 'FontWeight', 'bold', 'interpreter', 'latex');
ylabel('$F_X(x)$', 'FontName', 'Times New Roman', ...
	'FontSize', fontSize_, 'FontWeight', 'bold', 'interpreter', 'latex');
ylim([0 0.5])
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = fontSize_;

exportgraphics(ax, 'binomial2.png', 'Resolution', 300);



%%
clear variables; close all; clc;

fontSize_ = 40;

n_	= 10;
p_	= 0.05;
lambda_		= n_*p_;
for k_ = 1:n_
	FXBin(k_)	= nchoosek(n_, k_)* p_^k_ * (1 - p_)^(n_ - k_);
	FXPoisson(k_)	= (lambda_^k_ / factorial(k_) ) * exp(-lambda_);
end




figure('units', 'normalized', 'OuterPosition', [0.05 0.05 0.8 0.9]);
bar(1:n_, FXPoisson, 'EdgeColor','white', 'FaceColor',[0.8500 0.3250 0.0980], ...
	'FaceAlpha', 1, 'BarWidth', 0.6); hold on;
bar(1:n_, FXBin, 'EdgeColor','white', 'FaceColor','black', ...
	'FaceAlpha', 0.6, 'BarWidth', 0.3); hold on;

xlabel('Outcome $x$', 'FontName', 'Times New Roman', ...
	'FontSize', fontSize_, 'FontWeight', 'bold', 'interpreter', 'latex');
ylabel('$F_X(x)$', 'FontName', 'Times New Roman', ...
	'FontSize', fontSize_, 'FontWeight', 'bold', 'interpreter', 'latex');
ylim([0 0.5])
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = fontSize_;

exportgraphics(ax, 'binomial-poisson.png', 'Resolution', 300);