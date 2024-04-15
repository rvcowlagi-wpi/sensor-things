clear variables; close all; clc

positionX0	= [10; 10];		% km
timeStamps	= 0:0.01:1;		% hr
nTimeStamps = length(timeStamps);

velocityV	= [50; 0];		% km/hr
positionX(1, :)	= positionX0(1) + timeStamps*velocityV(1);
positionX(2, :)	= positionX0(2) + timeStamps*velocityV(2);

rangeTrue	= (positionX(1, :).^2 + positionX(2, :).^2).^(0.5);
bearTrue	= atan2(positionX(2, :), positionX(1, :));

vNoiseRange	= (0.3)*randn(1, length(timeStamps));							% Range error std. dev. is 0.3 km
vNoiseBear	= (2*pi/180)*randn(1, length(timeStamps));						% Bearing error std. dev. is 2 deg
zRange		= rangeTrue + vNoiseRange;
zBear		= bearTrue + vNoiseBear;

vRangeClutter	= 0.5*randn(1, length(timeStamps));
vBearClutter= (2*pi/180)*randn(1, length(timeStamps));						% Bearing error std. dev. is 2 deg
zRange2		= rangeTrue + vRangeClutter;
zBear2		= bearTrue + vBearClutter;


vRangeClutter	= 1*randn(1, length(timeStamps));
vBearClutter= (10*pi/180)*randn(1, length(timeStamps));						% Bearing error std. dev. is 2 deg
zRange3		= rangeTrue + vRangeClutter;
zBear3		= bearTrue + vBearClutter;


brng0		= atan2(positionX0(2), positionX0(1));
range0		= norm(positionX0);
xi0			= [range0; velocityV(1)*cos(brng0); brng0; 0];
[~, xiSim]	= ode45(@(t,x) polar_kinematics_2D(t,x, velocityV(1)), timeStamps, xi0);


save data_2D_RT_tracker_clutter.mat zRange zRange2 zRange3 ...
	zBear zBear2 zBear3 timeStamps rangeTrue bearTrue

% fig1= figure;
% plot(timeStamps, rangeTrue, 'LineWidth', 3);% , timeStamps, zRange);
% hold on; plot(timeStamps, xiSim(:, 1), 'LineWidth', 3)
% ax1	= gca;
% make_nice_figures(fig1, ax1, 18, [], 'Time', 'Range', 'Range', [],[],[],[]);
% 
% 
% fig2= figure;
% plot(timeStamps, bearTrue*180/pi, timeStamps, zBear*180/pi); hold on;
% plot(timeStamps, xiSim(:, 3)*180/pi, 'LineWidth', 3)
% ax2	= gca;
% make_nice_figures(fig2, ax2, 18, [], 'Time', 'Bearing (deg)', 'Bearing', [],[],[],[]);
% 
% figure;
% % plot(timeStamps, rangeTrue, 'LineWidth', 3);% , timeStamps, zRange);
% hold on; plot(timeStamps, xiSim(:, 3), 'LineWidth', 3)
% ax1	= gca;
% make_nice_figures(gcf, gca, 18, [], 'Time', 'Range Rate', 'Range Rate', [],[],[],[]);

figure;
for m1 = 1:nTimeStamps
	cla;
	plot(rangeTrue(m1)*cos(bearTrue(m1)), rangeTrue(m1)*sin(bearTrue(m1)), ...
		'.', 'MarkerSize', 20); hold on;
	ax = gca;
	plot(zRange(m1)*cos(zBear(m1)), zRange(m1)*sin(zBear(m1)), ...
		'x', 'MarkerSize', 10); ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
	plot(zRange2(m1)*cos(zBear2(m1)), zRange2(m1)*sin(zBear2(m1)), ...
		'x', 'MarkerSize', 10); ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
	plot(zRange3(m1)*cos(zBear3(m1)), zRange3(m1)*sin(zBear3(m1)), ...
		'x', 'MarkerSize', 10)
	make_nice_figures(gcf, gca, 18, [], '$y_1$', '$y_2$', ...
		'Range Rate', [],[],[0 80],[0 20]);
	drawnow()
	pause(0.1)
	
end

function xDot = polar_kinematics_2D(t_, x_, V_)
	xDot	= zeros(4, 1);

	xDot(1) = x_(2);
	xDot(3)	= x_(4);
% 	xDot(1)	= V_*cos(x_(3));
% 	xDot(3)	= -V_*sin(x_(3)) / x_(1);
	xDot(2)	= -V_*x_(4)*sin(x_(3));
	xDot(4)	= V_*x_(2)*sin(x_(3)) / (x_(1)^2) - V_*x_(4)*cos(x_(3)) / x_(1);
end