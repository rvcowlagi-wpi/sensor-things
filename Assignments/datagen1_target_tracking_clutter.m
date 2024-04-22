%% Asynchronous, multi-rate position, velocity, and time
clear variables; close all; clc


p0	= [8; 9];					% position, km
v0	= 40;						% speed, km/hr
psi0= 45*pi/180;				% heading angle, rad

dt_			= 0.001;
timeStamps	= 0:dt_:0.75;			% hr
nTimeStamps	= length(timeStamps);

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

acc			= 4*(10*sin(12*timeStamps) + 6*cos(7*timeStamps) - 2*sin(35*timeStamps));
steerRate	= 80*sin(35*timeStamps) + 50*cos(20*timeStamps);


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
groundTruthTimeStamps		= timeStamps;

groundTruthRange			= ( (groundTruthPosition(1, :)).^2 + ...
	(groundTruthPosition(2, :)).^2 ).^(0.5);
groundTruthBearing			= atan2( groundTruthPosition(2, :), groundTruthPosition(1, :) )*180/pi;

measuredRange1	= groundTruthRange(timeIndexRT) + stdDevRange*randn(1, nTimeStampsRT);
measuredBearing1	= groundTruthBearing(timeIndexRT) + stdDevBrng*randn(1, nTimeStampsRT);

measuredRange2	= groundTruthRange(timeIndexRT) + 2*stdDevRange*randn(1, nTimeStampsRT);
measuredBearing2	= groundTruthBearing(timeIndexRT) + 1.1*stdDevBrng*randn(1, nTimeStampsRT);

measuredRange3	= groundTruthRange(timeIndexRT) + 1.5*stdDevRange*randn(1, nTimeStampsRT);
measuredBearing3	= groundTruthBearing(timeIndexRT) + 1.5*stdDevBrng*randn(1, nTimeStampsRT);



% save assignment3_problem4.mat measuredRange measuredBearing ...
% 	stdDevRange stdDevBrng timeStampsRT
	
% save assignment4_problem2.mat measuredRange1 measuredBearing1 ...
% 	stdDevRange stdDevBrng timeStampsRT ...
% 	groundTruthTimeStamps groundTruthAcceleration ...
% 	groundTruthSpeed groundTruthHeading ... 
% 	groundTruthPosition groundTruthSteerRate

save assignment5_problem1.mat measuredRange1 measuredBearing1 ...
	measuredRange2 measuredBearing2 ...
	measuredRange3 measuredBearing3 ...
	stdDevRange stdDevBrng timeStampsRT ...
	groundTruthTimeStamps groundTruthAcceleration ...
	groundTruthSpeed groundTruthHeading ... 
	groundTruthPosition groundTruthSteerRate



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
function xDot = kinematics2D(t_, x_, u_, nTimeStamps_, dt)
		tIndex	= max(1, min(round(t_/dt), nTimeStamps_));
		ut_		= u_(:, tIndex);

		xDot	= [x_(3)*cos(x_(4)); x_(3)*sin(x_(4)); ut_(1); ut_(2)];
end