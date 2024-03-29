%% Synchronized position and time
clear variables; close all; clc

p0	= [10; 10];		% position, m
v0	= [1; 0];		% velocity, m/s

dt_			= 0.01;
timeStamps	= 0:0.01:10;		% s
nTimeStamps	= length(timeStamps);
acc			= zeros(2, nTimeStamps);
acc(1, :)	= sin(timeStamps) + 2*cos(2*timeStamps);
acc(2, :)	= 0.5 + 1*sin(4*timeStamps) + 1*cos(2.5*timeStamps);

[~, xSim]	= ode45(@(t,x) kinematics2D(t, x, acc, nTimeStamps, dt_), ...
	timeStamps, [p0; v0]);	

stdDevPos	= 1;
stdDevVel	= 0.5;
stdDevAcc	= 0.5;

measNoiseCovarPosn	= [stdDevPos 0; 0 stdDevPos];
measNoiseCovarAcc	= [stdDevAcc 0; 0 stdDevAcc];

RPosn	= chol(measNoiseCovarPosn);
RAcc	= chol(measNoiseCovarAcc);

groundTruthPosition			= xSim(:, 1:2)';
groundTruthVelocity			= xSim(:, 3:4)';
groundTruthAcceleration		= acc;

measuredPosition			= (xSim(:, 1:2))' + RPosn*randn(2, nTimeStamps);
measuredAcceleration		= acc + RAcc*randn(2, nTimeStamps);

save assignment3_problem1.mat measuredAcceleration measuredPosition ...
	groundTruthAcceleration groundTruthVelocity groundTruthPosition ...
	timeStamps nTimeStamps stdDevPos stdDevPos


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






%% Asynchronous, multi-rate position, velocity, and time
clear variables; close all; clc


p0	= [8; 3];		% position, m
v0	= [0.5; 0.8];		% velocity, m/s

dt_			= 0.001;
timeStamps	= 0:dt_:10;		% s
nTimeStamps	= length(timeStamps);
acc			= zeros(2, nTimeStamps);
acc(1, :)	= 0.7*(2*sin(5*timeStamps) + 4*cos(1.8*timeStamps));
acc(2, :)	= 0.6*sin(4*timeStamps) + 2.8*cos(2.5*timeStamps) - 0.3*cos(sin(5*timeStamps));

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

save assignment3_problem2.mat measuredAcceleration measuredPosition ...
	measuredVelocity stdDevPos stdDevAcc stdDevVel ...
	timeStampsAcc timeStampsVel timeStampsPosn
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



%% Synchronized position and time, no covariances

clear variables; close all; clc

p0	= [10; 10];		% position, m
v0	= [-1.5; 2];		% velocity, m/s

dt_			= 0.01;
timeStamps	= 0:dt_:10;		% s
nTimeStamps	= length(timeStamps);
acc			= zeros(2, nTimeStamps);
acc(1, :)	= 1.7*(2*sin(5*timeStamps) + 2.8*cos(1.8*timeStamps));
acc(2, :)	= 1.9*sin(5*timeStamps) + 3.1*cos(2.5*timeStamps) - 1*cos(sin(5*timeStamps));

[~, xSim]	= ode45(@(t,x) kinematics2D(t, x, acc, nTimeStamps, dt_), ...
	timeStamps, [p0; v0]);	

stdDevPos	= 2;
stdDevVel	= 0.5;
stdDevAcc	= 0.1;

measNoiseCovarPosn	= [stdDevPos 0; 0 stdDevPos];
measNoiseCovarAcc	= [stdDevAcc 0; 0 stdDevAcc];

RPosn	= chol(measNoiseCovarPosn);
RAcc	= chol(measNoiseCovarAcc);

groundTruthPosition			= xSim(:, 1:2)';
groundTruthVelocity			= xSim(:, 3:4)';
groundTruthAcceleration		= acc;

measuredPosition			= (xSim(:, 1:2))' + RPosn*randn(2, nTimeStamps);
measuredAcceleration		= acc + RAcc*randn(2, nTimeStamps);

save assignment3_problem3.mat measuredAcceleration measuredPosition ...
	timeStamps nTimeStamps


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