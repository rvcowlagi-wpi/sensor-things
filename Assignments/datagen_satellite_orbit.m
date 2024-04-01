function datagen_satellite_orbit()

%% Constants and parameters
close all; clc

dt_			= 60;						% s
timeStamps	= 0:dt_:(3*60*60);			% s
nTimeStamps	= length(timeStamps);

REarth		= 6378;				% km
MU_			= 398600;			% km^3 / s^2
H0_			= 65750;			% km^2 / s

% thtaPlot = 0:0.01:2*pi;
% plot(thtaPlot, (H0_^2 / MU_)*(1 + 0.6*cos(thtaPlot)), 'LineWidth', 2 )

% return

% Based on Ex 2.5 in Curtis
eccKepler	= 0.6;				% eccentricity
thta0		= 0;
r0			= (H0_^2 / MU_) / (1 + eccKepler*cos(thta0));
spd0		= H0_ / r0;			% km / s
thtaDot0	= spd0 / r0;		% rad/s
rDot0		= 0;

stddevProcNoiseRange	= sqrt(1E-8);
stddevProcNoiseBearing	= sqrt(1E-8);

vRange		= stddevProcNoiseRange*randn(1, nTimeStamps);
vThta		= stddevProcNoiseBearing*randn(1, nTimeStamps);

[~, xSim]	= ode45(@orbit_kinematics_planar, ...
	timeStamps, [r0; rDot0; thta0; thtaDot0]);

stdDevRange		= 30;			% km
stdDevRangeRate	= 1*pi/180;		% km/s

groundTruthRange			= xSim(:, 1)';
groundTruthRangeRate		= xSim(:, 2)';
groundTruthBearing			= xSim(:, 3)';
groundTruthBearingRate		= xSim(:, 4)';

timeIndexMeasured			= 1:1:nTimeStamps;
timeStampsMeasured			= timeStamps(timeIndexMeasured);
nTimeStampsMeasured			= length(timeStampsMeasured);

measuredRange		= groundTruthRange(timeIndexMeasured) + stdDevRange*randn(1, nTimeStampsMeasured);
measuredRangeRate	= groundTruthRangeRate(timeIndexMeasured) + stdDevRangeRate*randn(1, nTimeStampsMeasured);

save assignment4_problem3.mat measuredRange measuredRangeRate ...
	stdDevRange stdDevRangeRate timeStamps nTimeStamps ...
	groundTruthRange groundTruthRangeRate groundTruthBearing ...
	groundTruthBearingRate ...
	stddevProcNoiseRange stddevProcNoiseBearing


figure;
plot( xSim(:, 1).*cos(xSim(:, 3)), xSim(:, 1).*sin(xSim(:, 3)), 'LineWidth', 2 );
hold on; axis equal
% plot( measuredRange.*cos(measuredBearing), measuredRange.*sin(measuredBearing), '.'); 
make_nice_figures(gcf, gca, 14, 'Orbit', '$y_1$ (km)', '$y_2$ (km)', [], [], [], [], [])

%%
	function xDot = orbit_kinematics_planar(t_, x_)							% Process noise added in discrete model, not here
		tIndex	= max(1, min(round(t_/dt_), nTimeStamps));

		r_		= x_(1);
		rDot	= x_(2);
		thta_	= x_(3);
		thtaDot	= x_(4);

		vrt_	= vRange(tIndex);
		vtt_	= vThta(tIndex);

		xDot	= [...
			rDot; ...
			r_*(thtaDot^2) - MU_ / r_^2 + vrt_; ...
			thtaDot; ...
			-2*rDot*thtaDot/r_ + (1 / r_)*vtt_];
	end
end