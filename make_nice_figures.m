function make_nice_figures(fig_, ax_, fontSize_, xLabel_, yLabel_, figTitle_, figPosition_)

fig_.Units		= 'normalized';
f
fig_.Position	= [0.1 0.1 0.5*[1 1]];

figure('Units','normalized', ) 
subplot(211); hold on; plot(time_stamps, storeXHat(1,:), 'LineWidth', 2)
title('State estimate mean $\hat{x}_1$ (position)', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', fontSize_);
ylabel('$\hat{x}_1$ (m)', 'Interpreter', 'latex', 'FontSize', fontSize_);
ax = gca; ax.FontSize = fontSize_; ax.FontName = 'Times New Roman';