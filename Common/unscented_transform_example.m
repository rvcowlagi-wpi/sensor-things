clear variables; close all; clc


yTrue		= [0; 1];
stdDevRange = 0.02;
stdDevBrng	= 15*pi/180;

rangeTrue	= norm(yTrue);
brngTrue	= atan2(yTrue(2), yTrue(1));

nSamples	= 1E3;

zRange		= rangeTrue + stdDevRange*randn(nSamples, 1);
zBrng		= brngTrue + stdDevBrng*randn(nSamples, 1);

y1PlotSamples	= zRange.*cos(zBrng);
y2PlotSamples	= zRange.*sin(zBrng);

figure;
plot(y1PlotSamples, y2PlotSamples, 'x'); axis equal; hold on;
make_nice_figures(gcf, gca, 18, [], '$y_1$', '$y_2$', [], [], [], [-1, 1], [0.8 1.2])

sampleMean	= mean([y1PlotSamples y2PlotSamples]);
sampleVar1	= var(y1PlotSamples);
sampleVar2	= var(y2PlotSamples);

yBar		= sampleMean;
p11			= sampleVar1;
p22			= sampleVar2;

disp([sqrt(p11) sqrt(p22)])





y1PlotTmp	= linspace( (yBar(1) - sqrt(p11)), (yBar(1) + sqrt(p11)), 100);
y2PlotTmp1	= sqrt(p22)*(1 - ((y1PlotTmp - yBar(1))/sqrt(p11)).^2).^0.5 + yBar(2);
y2PlotTmp2	= -sqrt(p22)*(1 - ((y1PlotTmp - yBar(1))/sqrt(p11)).^2).^0.5 + yBar(2);
plot(yBar(1), yBar(2), '.', 'MarkerSize', 20);
ax = gca; ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot(y1PlotTmp, y2PlotTmp1, 'LineWidth', 3);
ax = gca; ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot(y1PlotTmp, y2PlotTmp2, 'LineWidth', 3)

rangeMean	= mean(zRange);		% Sanity check
brngMean	= mean(zBrng);

P0			= [stdDevRange^2 0; 0 stdDevBrng^2];
jacA		= [...
	cos(brngMean) -rangeMean*sin(brngMean); ...
	sin(brngMean) -rangeMean*cos(brngMean)];
P1			= jacA * P0 * jacA'

yBar		= rangeMean*[cos(brngMean); sin(brngMean)];
p11			= P1(1,1);
p22			= P1(2,2);

disp([sqrt(p11) sqrt(p22)])


y1PlotTmp	= linspace( (yBar(1) - sqrt(p11)), (yBar(1) + sqrt(p11)), 100);
y2PlotTmp1	= sqrt(p22)*(1 - ((y1PlotTmp - yBar(1))/sqrt(p11)).^2).^0.5 + yBar(2);
y2PlotTmp2	= -sqrt(p22)*(1 - ((y1PlotTmp - yBar(1))/sqrt(p11)).^2).^0.5 + yBar(2);
plot(yBar(1), yBar(2), '.', 'MarkerSize', 20);
ax = gca; ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot(y1PlotTmp, y2PlotTmp1, 'LineWidth', 3);
ax = gca; ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot(y1PlotTmp, y2PlotTmp2, 'LineWidth', 3)


%%
nX	= 2;
S	= chol(nX*P0);

xSigma(:, 1)	= [rangeMean; brngMean] + S(1, :)';
xSigma(:, 2)	= [rangeMean; brngMean] + S(2, :)';
xSigma(:, 3)	= [rangeMean; brngMean] - S(1, :)';
xSigma(:, 4)	= [rangeMean; brngMean] - S(2, :)';

