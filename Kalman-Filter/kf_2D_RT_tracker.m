%% Kalman Filter implementation for 2D Range-Bearing Tracking
% Copy-pasted from 1D navigation; fix comments later

clear variables; close all; clc

%%
% Load given data file; it must be in the same folder as this file.
% In practice, data are not available prerecorded like this example.
% Instead, data are received in real time.
% This data file loads the following variables with self-evident meanings:
% |timeStamps|, |accelerometer|, |gps_position|, and |gps_speed.|
load data_2D_RT_tracker.mat

nStates = 4;	% States are r, rDot, theta, thetaDot
nMeas	= 2;	% Measurements are r, theta

%% Error Covariances
Q	= [(100)^2 0; 0 (20*pi/180)^2];
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
% Remember that the predictive model is 
%%
% $$ \underline{x}(k) = F\underline{x}(k-1) + G_1\underline{u}(k-1) + \underline{w}(k-1)$$
%%
% In this example $\underline{u}$ and $\underline{w}$ are scalars.
A	= [ [1 dt; 0 1] zeros(2) ; zeros(2) [1 dt; 0 1] ];
B	= zeros(4,1);
G	= [ [0; dt] zeros(2, 1); zeros(2, 1) [0; dt] ];
%%
% The measurement model is $\underline{z} = C\underline{x}.$
C	= [1 0 0 0; 0 0 1 0];

%% Initialization
%%
% Set initial values for state estimate $\underline{\hat{x}}$ and
% estimation error covariance $P$.
xHat	= [zRange(1); 0; zBear(1); 0];	
P		= diag([R(1,1), 0.1, R(2,2), 0]);

[PSteadyState,~,~]	= idare(A', C', (G*Q*G'), R, [], []);
PSteadyStateStore	= reshape(PSteadyState, nStates^2, 1);
LSteadyState		= PSteadyState * C' / (C*PSteadyState*C' + R);
%% 
% Memory allocation for recording the state estimate and error covariance
% over time. We will store the result of each iteration in the variables
% |xhat_store| and |P_store|. We also store the trace of the P matrix at
% each iteration in |P_trace_store|. Each column of these variables stores
% $\underline{\hat{x}}$ and $P$ at each iteration. 
storeXHat	= zeros(nStates, nTimeStamps);
storeP		= zeros(nStates^2, nTimeStamps);
storePTrace	= zeros(1, nTimeStamps);
storeInnov	= zeros(nMeas, nTimeStamps);
%%
% Therefore, the first column stores the initial values. In general the
% $k^\mathrm{th}$ column of |xhat_store| stores $\underline{\hat{x}}(k-1)$
% and so on. The $(k-1)$ is because MATLAB's array numbering starts at 1,
% whereas our $k$ values (time stamps) start at 0.
storeXHat(:, 1)		= xHat;
storePTrace(:, 1)	= trace(P);
storeP(:, 1)		= reshape(P, nStates^2, 1); 
%%
% The |reshape| command in the previous line stores the matrix P as a 4x1
% array for convenience. Otherwise we would need a three-dimensional array.

%% Run Kalman Filter Iterations
%%
% There are as many iterations as time stamps, minus 1 because we
% initialize using data at the first time stamp data. Think of the
% iteration variable |m1| = $k-1.$
for m1 = 1:(nTimeStamps-1)
%%
% Prediction equations to get a preliminary estimate and error covariance.
	u			= 0;
	xHatMinus	= A*xHat + B*u;
	PMinus		= A*P*A' + G*Q*G';
%%	
% Compute Kalman gain. Note the use of |/| command instead of inverse.
	L			= (PMinus * C') / (C * PMinus * C' + R );
	
%% 
% Correction equations to get new estimate and error covariance at this
% iteration. Note the |eye(2)| command, which is a 2x2 identity matrix.
	z			= [zRange(m1 + 1); zBear(m1 + 1)];
	thisInnov	= z - C*xHatMinus;
	xHat		= xHatMinus + L*(thisInnov);
% 	xhat		= x_minus + LSteadyState*thisInnovation;
	P			= (eye(nStates) - L*C)*PMinus;
	
%% 
% Store the newly computed estimate and error covariance.
	storeXHat(:, m1+1)		= xHat;
	storeP(:, m1+1)			= reshape(P, nStates^2, 1);
	storePTrace(:, m1+1)	= trace(P);
	storeInnov(:, m1+1)		= thisInnov;
end
%%
% End of Kalman filter iterations.

%% Innovation autocorrelation
innovCovar		= C * PSteadyState * C' + R;

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
plot(timeStamps, 2*sqrt(innovCovar(1,1))*ones(1, length(timeStamps)), ...
	'LineWidth', 2, 'LineStyle', '--'); ax.ColorOrderIndex = 1;
plot(timeStamps, -2*sqrt(innovCovar(1,1))*ones(1, length(timeStamps)), ...
	'LineWidth', 2, 'LineStyle', '--')

subplot(122)
plot(timeStamps, storeInnov(2, :), 'LineWidth', 2); hold on;
ax = gca; ax.ColorOrderIndex = 2;
plot(timeStamps, 2*sqrt(innovCovar(2,2))*ones(1, length(timeStamps)), ...
	'LineWidth', 2, 'LineStyle', '--'); ax.ColorOrderIndex = 2;
plot(timeStamps, -2*sqrt(innovCovar(2,2))*ones(1, length(timeStamps)), ...
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
