%% Asynchronous, multi-rate position, velocity, and time
clear variables; close all; clc


p0	= [8; 3];			% position, km
v0	= [0; 18];		% velocity, km/hr

dt_			= 0.001;
timeStamps	= 0:dt_:1;		% hr
nTimeStamps	= length(timeStamps);
acc			= zeros(2, nTimeStamps);
for m1 = 100:100:nTimeStamps
	acc(:, m1 - 99 : m1) = kron(40*randn(2, 1), ones(1, 100));
end


[~, xSim]	= ode45(@(t,x) kinematics2D(t, x, acc, nTimeStamps, dt_), ...
	timeStamps, [p0; v0]);	

stdDevPos	= 1;
stdDevVel	= 0.5;
stdDevAcc	= 0.5;

timeStampsPosn	= timeStamps(2:100:end);
timeStampsVel	= timeStamps(3:100:end);
timeStampsAcc	= timeStamps(1:10:end);


measNoiseCovarPosn	= [stdDevPos 0; 0 stdDevPos];
measNoiseCovarVel	= [stdDevVel 0; 0 stdDevVel];
measNoiseCovarAcc	= [stdDevAcc 0; 0 stdDevAcc];

RPosn	= chol(measNoiseCovarPosn);
RVel	= chol(measNoiseCovarVel);
RAcc	= chol(measNoiseCovarAcc);


groundTruthPosition			= xSim(:, 1:2)';
groundTruthVelocity			= xSim(:, 3:4)';
groundTruthAcceleration		= acc;

measuredPosition			= (xSim(2:100:end, 1:2))' + RPosn*randn(2, length(timeStampsPosn));
measuredVelocity			= (xSim(3:100:end, 3:4))' + RVel*randn(2, length(timeStampsVel));
measuredAcceleration		= acc(:, 1:10:end) + RAcc*randn(2, length(timeStampsAcc));

% save assignment3_problem2.mat measuredAcceleration measuredPosition ...
% 	measuredVelocity stdDevPos stdDevAcc stdDevVel ...
% 	timeStampsAcc timeStampsVel timeStampsPosn
% 	groundTruthAcceleration groundTruthVelocity groundTruthPosition ...
	


figure;
subplot(311); plot(timeStamps, xSim(:, 1), timeStamps, xSim(:, 2), 'LineWidth', 2);
make_nice_figures(gcf, gca, 14, 'Position', ...
	'Time (s)', '$\mathbf{p}$ (m)', [], [], [], [], [])

subplot(312); plot(timeStamps, xSim(:, 3), timeStamps, xSim(:, 4), 'LineWidth', 2);
make_nice_figures(gcf, gca, 14, 'Velocity', ...
	'Time (s)', '$\mathbf{\dot{p}}$ (m/s)', [], [], [], [], [])

subplot(313); plot(timeStamps, acc(1, :), timeStamps, acc(2, :), 'LineWidth', 2);
make_nice_figures(gcf, gca, 14, 'Acceleration', ...
	'Time (s)', '$\mathbf{a}$ (m/s$^2$)', [], [], [], [], [])


figure;
plot(xSim(:, 1), xSim(:, 2), 'LineWidth', 2); axis equal
make_nice_figures(gcf, gca, 14, 'Position', '$x$ (m)', '$y$ (m)', [], [], [], [], [])



%%
function xDot = kinematics2D(t_, x_, u_, nTimeStamps_, dt)
		tIndex	= max(1, min(round(t_/dt), nTimeStamps_));
		ut_		= u_(:, tIndex);
		xDot	= [x_(3:4); ut_];
end