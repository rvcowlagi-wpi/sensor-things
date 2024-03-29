function datagen_sat()

%% Constants and parameters
close all; clc

dt_			= 0.1;
timeStamps	= 0:dt_:(0.01*60*60);			% s
nTimeStamps	= length(timeStamps);

REarth		= 6378;				% km
MU_			= 398600;			% km^3 / s^2

% Based on Ex 2.5 in Curtis
eccKepler	= 0.6;				% eccentricity
r0			= 400 + REarth;		% perigee, km
thta0		= 0;
spd0		= 9.7;				% km / s
h0			= 65750;			% km^2 / s
thtaDot0	= spd0 / r0		% rad/s
rDot0		= 0;

[~, xSim]	= ode45(@orbit_kinematics_planar, ...
	timeStamps, [r0; rDot0; thta0; thtaDot0]);

plot(timeStamps, xSim(:, 1), 'LineWidth', 2)

return


acc			= zeros(1, nTimeStamps);	% tangential acceleration, km/hr^2
steerRate	= zeros(1, nTimeStamps);	% steering rate, rad/hr

% accStep		= 100:100:nTimeStamps;
% for m1 = accStep
% 	acc(:, (m1 - accStep(1) + 1):m1) = ...
% 		kron(20*randn(1, 1), ones(1, accStep(1)));
% 
% 	steerRate(:, (m1 - accStep(1) + 1):m1) = ...
% 		kron(10*2*pi*randn(1, 1), ones(1, accStep(1)));
% end

acc			= 8*(10*sin(18*timeStamps) + 6*cos(15*timeStamps) - 2*sin(50*timeStamps));
steerRate	= 80*sin(80*timeStamps) + 50*cos(20*timeStamps);


[~, xSim]	= ode45(@(t,x) kinematics2D(t, x, [acc; steerRate], nTimeStamps, dt_), ...
	timeStamps, [p0; v0; psi0]);	

stdDevRange	= 0.3;			% km
stdDevBrng	= 5;			% deg

timeIndexRT		= 1:15:nTimeStamps;
timeStampsRT	= timeStamps(timeIndexRT);
nTimeStampsRT	= length(timeStampsRT);

groundTruthPosition			= xSim(:, 1:2)';
groundTruthSpeed			= xSim(:, 3)';
groundTruthHeading			= xSim(:, 4)';
groundTruthAcceleration		= acc;
groundTruthSteerRate		= steerRate;

groundTruthRange			= ( (groundTruthPosition(1, :)).^2 + ...
	(groundTruthPosition(2, :)).^2 ).^(0.5);
groundTruthBearing			= atan2( groundTruthPosition(2, :), groundTruthPosition(1, :) )*180/pi;

measuredRange	= groundTruthRange(timeIndexRT) + stdDevRange*randn(1, nTimeStampsRT);
measuredBearing	= groundTruthBearing(timeIndexRT) + stdDevBrng*randn(1, nTimeStampsRT);

save assignment4_problem2.mat measuredRange measuredBearing ...
	stdDevRange stdDevBrng timeStampsRT
	

figure;
subplot(211); plot(timeStamps, xSim(:, 1), timeStamps, xSim(:, 2), 'LineWidth', 2);
make_nice_figures(gcf, gca, 14, 'Position', ...
	'Time (s)', 'Position (km)', [], [], [], [], [])

subplot(212); plot(timeStamps, xSim(:, 3), 'LineWidth', 2);
make_nice_figures(gcf, gca, 14, 'Velocity', ...
	'Time (s)', 'Speed (km/hr)', [], [], [], [], [])

figure;
subplot(211); plot(timeStamps, acc, 'LineWidth', 2);
make_nice_figures(gcf, gca, 14, 'Acceleration', ...
	'Time (s)', 'Accel.', [], [], [], [], [])
subplot(212); plot(timeStamps, steerRate, 'LineWidth', 2);
make_nice_figures(gcf, gca, 14, 'Acceleration', ...
	'Time (s)', 'Steer rate', [], [], [], [], [])


figure;
plot(xSim(:, 1), xSim(:, 2), 'LineWidth', 2); axis equal
make_nice_figures(gcf, gca, 14, 'Position', '$x$ (km)', '$y$ (km)', [], [], [], [], [])



%%
	function xDot = orbit_kinematics_planar(t_, x_)							% Process noise added in discrete model, not here
% 		tIndex	= max(1, min(round(t_/dt), nTimeStamps));

		r_		= x_(1);
		rDot	= x_(2);
		thta_	= x_(3);
		thtaDot	= x_(4);

		xDot	= [r_; r_*(thtaDot^2) - MU_ / r_^2; thtaDot; -2*rDot*thtaDot/r_];
	end
end