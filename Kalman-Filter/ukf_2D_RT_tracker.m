%% Kalman Filter implementation for 2D Range-Bearing Tracking
% Copy-pasted from 1D navigation; fix comments later

clear variables; close all; clc

%%
load data_2D_RT_tracker.mat

nStates = 4;	% States are r, rDot, theta, thetaDot
nMeas	= 2;	% Measurements are r, theta

V		= 50;

%% Error Covariances
Q	= 20*[(1E-2)^2 0; 0 (1E-2*pi/180)^2];
R	= [0.3^2 0; 0 (2*pi/180)^2];

%% Time Step and Interval of Interest
%%
% We have uniformly spaced time stamps, so calculating the difference
% between any two successive instants gives us $dt.$
dt	= timeStamps(2) - timeStamps(1);
%%
% Note the number of time stamps, which will be the number of iterations
% for which the Kalman filter will run.
nTimeStamps = length(timeStamps);

%% Predictive Model and Measurement Model
G	= [ [0; dt] zeros(2, 1); zeros(2, 1) [0; dt] ];
C	= [1 0 0 0; 0 0 1 0];

%% Initialization
xHat	= [zRange(1); V; zBear(1); -V/zRange(1)];	
P		= diag([R(1,1), 0.1, R(2,2), 0]);
% 
% [PSteadyState,~,~]	= idare(A', C', (G*Q*G'), R, [], []);
% PSteadyStateStore	= reshape(PSteadyState, nStates^2, 1);
% LSteadyState		= PSteadyState * C' / (C*PSteadyState*C' + R);
%% 
storeXHat	= zeros(nStates, nTimeStamps);
storeP		= zeros(nStates^2, nTimeStamps);
storePTrace	= zeros(1, nTimeStamps);
storeInnov	= zeros(nMeas, nTimeStamps);
storeInnovCovar	= zeros(nMeas^2, nTimeStamps);

%%
storeXHat(:, 1)		= xHat;
storePTrace(:, 1)	= trace(P);
storeP(:, 1)		= reshape(P, nStates^2, 1); 

%% Run Kalman Filter Iterations
for m1 = 1:(nTimeStamps-1)

	%----- Prediction equations to get a preliminary estimate and error covariance.
	u			= 0;
	A			= eye(nStates) + jacobianA(xHat, V)*dt;

% 	xHatMinus	= xHat + polar_kinematics_2D([], xHat, V)*dt;				% First-order Euler approximation
	[~, xSim]	= ode45(@(t,x) polar_kinematics_2D(t,x, V), [0 dt], xHat);	% RK4
	xHatMinus	= xSim(end, :)';
	PMinus		= A*P*A' + G*Q*G';

	%----- Compute Kalman gain. Note the use of |/| command instead of inverse.
	L			= (PMinus * C') / (C * PMinus * C' + R );
	
	%----- Measurement update
	z			= [zRange(m1 + 1); zBear(m1 + 1)];
	thisInnov	= z - C*xHatMinus;
	
	xHat		= xHatMinus + L*(thisInnov);
	P			= (eye(nStates) - L*C)*PMinus;

	thisInnovCovar		= C * P * C' + R;
	
	%----- Store results
	storeXHat(:, m1+1)		= xHat;
	storeP(:, m1+1)			= reshape(P, nStates^2, 1);
	storePTrace(:, m1+1)	= trace(P);
	storeInnov(:, m1+1)		= thisInnov;
	storeInnovCovar(:, m1+1)= reshape(thisInnovCovar, nMeas^2, 1);
end

%% Innovation autocorrelation
nTimeStamps		= length(timeStamps);
nMeasDim		= 2;

innovAutoCorr	= zeros(1, nTimeStamps);
shiftedInnov	= zeros(nMeasDim, nTimeStamps);
twoSigmaBd		= zeros(1, nTimeStamps);
for m1 = 1:nTimeStamps - 10
	shiftedInnov(:,m1+1:end)= storeInnov(:, 2:(nTimeStamps - m1 + 1));
	innovAutoCorr(m1)		= (1 / (nTimeStamps - m1) ) * (...
			sum( shiftedInnov(1, m1+1:end).*storeInnov(1, m1+1:end) ) + ...
			sum( shiftedInnov(2, m1+1:end).*storeInnov(2, m1+1:end) ) );
	twoSigmaBd(m1)			= 2/sqrt(nTimeStamps - m1); 
end

%% Plot Results
fig1 = figure; 
subplot(221); hold on; plot(timeStamps, storeXHat(1,:), 'LineWidth', 2);
hold on;
plot(timeStamps, rangeTrue, 'LineWidth',  2)
make_nice_figures(gcf, gca, 18, [], 'Time (h)', ...
	'Range $\hat{x}_1 = r$ (km)', 'Range and Rate', [],[],[],[]);

subplot(222); hold on; plot(timeStamps, storeXHat(2,:), 'LineWidth', 2)
make_nice_figures(gcf, gca, 18, [], 'Time (h)', ...
	'Range rate $\hat{x}_2 = \dot{r}$ (km/hr)', 'Range and Rate', [],[],[],[]);

subplot(223); hold on; plot(timeStamps, storeXHat(3,:)*180/pi, 'LineWidth', 2);
hold on;
plot(timeStamps, bearTrue*180/pi, 'LineWidth',  2)
make_nice_figures(gcf, gca, 18, [], 'Time (h)', ...
	'Range $\hat{x}_3 = \theta$ (deg)', 'Bearing and Rate', [],[],[],[]);

subplot(224); hold on; plot(timeStamps, storeXHat(4,:)*180/pi, 'LineWidth', 2)
make_nice_figures(gcf, gca, 18, [], 'Time (h)', ...
	'Range rate $\hat{x}_4 = \dot{\theta}$ (deg/hr)', 'Bearing and Rate', [],[],[],[]);




fig2 = figure;
plot(timeStamps, storePTrace, 'LineWidth', 2);
make_nice_figures(gcf, gca, 18, [], 'Time (h)', 'tr$(P)$', 'Trace', [],[],[],[]);

fig3 = figure;
plot(timeStamps, storeP, 'LineWidth', 2); hold on;
make_nice_figures(gcf, gca, 18, [], 'Time (h)', '$p_{ij}$', 'Trace', [],[],[],[]);

ylabel('$p_{ij}$', 'Interpreter', 'latex', 'FontSize', 12);


figure;
subplot(121)
plot(timeStamps, storeInnov(1, :), 'LineWidth', 2); hold on;
ax = gca; ax.ColorOrderIndex = 1;
plot(timeStamps, 2*(storeInnovCovar(1,:).^0.5), ...
	'LineWidth', 2, 'LineStyle', '--'); ax.ColorOrderIndex = 1;
plot(timeStamps, -2*(storeInnovCovar(1,:).^0.5), ...
	'LineWidth', 2, 'LineStyle', '--')

subplot(122)
plot(timeStamps, storeInnov(2, :), 'LineWidth', 2); hold on;
ax = gca; ax.ColorOrderIndex = 2;
plot(timeStamps, 2*(storeInnovCovar(4, :).^0.5), ...
	'LineWidth', 2, 'LineStyle', '--'); ax.ColorOrderIndex = 2;
plot(timeStamps, 2*(storeInnovCovar(4, :).^0.5), ...
	'LineWidth', 2, 'LineStyle', '--')

title('Innovations sequence', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Innovation', 'Interpreter', 'latex', 'FontSize', 12);


figure;
plot(timeStamps, innovAutoCorr, 'LineWidth', 2); hold on;
plot(timeStamps, (2/sqrt(nTimeStamps))*ones(size(timeStamps)), 'LineWidth', 2, 'LineStyle', '--' );
ax = gca; ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot(timeStamps, (-2/sqrt(nTimeStamps))*ones(size(timeStamps)), 'LineWidth', 2, 'LineStyle', '--' )
title('Innovations autocorrelation', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Innovation', 'Interpreter', 'latex', 'FontSize', 12);


%%
function xDot = polar_kinematics_2D(t_, x_, V_)
	xDot	= zeros(4, 1);

% 	xDot(1) = x_(2);
% 	xDot(3)	= x_(4);
	xDot(1)	= V_*cos(x_(3));
	xDot(3)	= -V_*sin(x_(3)) / x_(1);
	xDot(2)	= -V_*x_(4)*sin(x_(3));
	xDot(4)	= V_*x_(2)*sin(x_(3)) / (x_(1)^2) - V_*x_(4)*cos(x_(3)) / x_(1);
end

function xSigma = generate_sigma_points(xHat_, PX_)
	lam0 = 1/3;


end