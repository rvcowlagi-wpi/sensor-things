clear variables; close all; clc

nSamples = 1E4;
y_		= zeros(nSamples, 1);
x_		= zeros(nSamples, 1);
xPlot	= 0.1:0.01:4;
for m1 = 1:nSamples
	x_(m1) = randn;
end

figure
histogram(x_, 'Normalization', 'pdf'); hold on;
plot(xPlot, normpdf(xPlot), 'LineWidth', 3)


figure;
y_		= x_.^2;
histogram(y_, 'Normalization','pdf'); hold on;
plot(xPlot, chi2pdf(xPlot, 1), 'LineWidth', 3)