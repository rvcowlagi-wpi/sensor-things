%% Kalman Filter implementation for 2D Range-Bearing Tracking
% Copy-pasted from 1D navigation; fix comments later

clear variables; close all; clc

%%
% Load given data file; it must be in the same folder as this file.
% In practice, data are not available prerecorded like this example.
% Instead, data are received in real time.
% This data file loads the following variables with self-evident meanings:
% |time_stamps|, |accelerometer|, |gps_position|, and |gps_speed.|
load data_2D_RT_tracker.mat

%% Error Covariances
Q	= [(0.5)^2 0; 0 (5*pi/180)^2];
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
A	= [1 dt; 0 1];
B	= [0; 0];
G	= [0; dt];
%%
% The measurement model is $\underline{z} = C\underline{x}.$
C	= [1 0];

%% Initialization
%%
% Set initial values for state estimate $\underline{\hat{x}}$ and
% estimation error covariance $P$.
xHat	= [rangeMeas(1); 0];
P		= [R 0; 0 1];

[PSteadyState,~,~]	= idare(A', C', (G*Q*G'), R, [], []);
PSteadyStateStore	= reshape(PSteadyState, 4, 1);
LSteadyState		= PSteadyState * C' / (C*PSteadyState*C' + R);
%% 
% Memory allocation for recording the state estimate and error covariance
% over time. We will store the result of each iteration in the variables
% |xhat_store| and |P_store|. We also store the trace of the P matrix at
% each iteration in |P_trace_store|. Each column of these variables stores
% $\underline{\hat{x}}$ and $P$ at each iteration. 
storeXHat	= zeros(2, nTimeStamps);
storeP		= zeros(4, nTimeStamps);
storePTrace	= zeros(1, nTimeStamps);
storeInnov	= zeros(2, nTimeStamps);
%%
% Therefore, the first column stores the initial values. In general the
% $k^\mathrm{th}$ column of |xhat_store| stores $\underline{\hat{x}}(k-1)$
% and so on. The $(k-1)$ is because MATLAB's array numbering starts at 1,
% whereas our $k$ values (time stamps) start at 0.
storeXHat(:, 1)		= xHat;
storePTrace(:, 1)	= trace(P);
storeP(:, 1)		= reshape(P, 4, 1); 
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
	z			= rangeMeas(m1 + 1);
	thisInnov	= z - C*xHatMinus;
	xHat		= xHatMinus + L*(thisInnov);
% 	xhat		= x_minus + LSteadyState*thisInnovation;
	P			= (eye(2) - L*C)*PMinus;
	
%% 
% Store the newly computed estimate and error covariance.
	storeXHat(:, m1+1)		= xHat;
	storeP(:, m1+1)			= reshape(P, 4, 1);
	storePTrace(:, m1+1)	= trace(P);
	storeInnov(:, m1+1)		= thisInnov;
end
%%
% End of Kalman filter iterations.

%% Innovation autocorrelation
nMeasDim		= 2;

innovAutoCorr	= zeros(1, nTimeStamps);
shiftedInnov	= zeros(nMeasDim, nTimeStamps);
for m1 = 1:(nTimeStamps - 75)
	shiftedInnov(:,m1+1:end)= storeInnov(:, 2:(nTimeStamps - m1 + 1));
	innovAutoCorr(m1)		= (1 / (nTimeStamps - m1) ) * (...
			sum( shiftedInnov(1, m1+1:end).*storeInnov(1, m1+1:end) ) + ...
			sum( shiftedInnov(2, m1+1:end).*storeInnov(2, m1+1:end) ) );
end

%% Plot Results
%% 
% Note that every plot has a title and has both axes labeled with units.
% This is standard practice. Plots are meaningless without labeled axes.
% All plots in your assignments should be similarly labeled. Units for the
% trace of $P$ are tricky because we did not non-dimensionalize the problem
% (as is often done in practice). Therefore, we are adding numbers in
% different units: m and m/s, which is not quite correct, but suffices to
% illustrate in this example that the filter works as expected. The
% decreasing trace of $P$ indicates increasing confidence in the estimate.
% We cannot achieve $tr(P) = 0$ because sensor noise is always present.
figure; 
subplot(211); hold on; plot(timeStamps, storeXHat(1,:), 'LineWidth', 2)
title('State estimate mean $\hat{x}_1$ (position)', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\hat{x}_1$ (m)', 'Interpreter', 'latex', 'FontSize', 12);

subplot(212); hold on; plot(timeStamps, storeXHat(2,:), 'LineWidth', 2)
title('State estimate mean $\hat{x}_2$ (speed)', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\hat{x}_2$ (m/s)', 'Interpreter', 'latex', 'FontSize', 12);


figure;
plot(timeStamps, storePTrace, 'LineWidth', 2)
title('Trace of estimation error covariance $P$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('tr$(P)$, units N/A', 'Interpreter', 'latex', 'FontSize', 12);


figure;
plot(timeStamps, storeP, 'LineWidth', 2); hold on;
ax = gca;
ax.ColorOrderIndex = 1;
plot(timeStamps, kron(PSteadyStateStore, ones(1, length(timeStamps))), 'LineWidth', 2, 'LineStyle', '--')
title('Elements of estimation error covariance $P$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$p_{ij}$', 'Interpreter', 'latex', 'FontSize', 12);
legend('$p_{11}$', '$p_{12}$', '$p_{21}$', '$p_{22}$', 'Interpreter', 'latex' )


figure;
plot(timeStamps, storeInnov, 'LineWidth', 2)
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
