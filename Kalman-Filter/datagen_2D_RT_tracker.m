clear variables; close all; clc

positionX0	= [10; 10];		% km
timeStamps	= 0:0.01:1;		% hr

velocityV	= [300; 0];		% km/hr
positionX(1, :)	= positionX0(1) + timeStamps*velocityV(1);
positionX(2, :)	= positionX0(2) + timeStamps*velocityV(2);

rangeTrue	= (positionX(1, :).^2 + positionX(2, :).^2).^(0.5);
bearTrue	= atan2(positionX(2, :), positionX(1, :));

vNoiseRange	= (0.3)*randn(1, length(timeStamps));							% Range error std. dev. is 0.3 km
vNoiseBear	= (2*pi/180)*randn(1, length(timeStamps));						% Bearing error std. dev. is 2 deg
zRange		= rangeTrue + vNoiseRange;
zBear		= bearTrue + vNoiseBear;


save data_2D_RT_tracker.mat zRange zBear timeStamps rangeTrue bearTrue

fig1= figure;
plot(timeStamps, rangeTrue, timeStamps, zRange);
ax1	= gca;
make_nice_figures(fig1, ax1, 18, [], 'Time', 'Range', 'Range', [],[],[],[]);


fig2= figure;
plot(timeStamps, bearTrue*180/pi, timeStamps, zBear*180/pi);
ax2	= gca;
make_nice_figures(fig2, ax2, 18, [], 'Time', 'Bearing (deg)', 'Bearing', [],[],[],[]);
