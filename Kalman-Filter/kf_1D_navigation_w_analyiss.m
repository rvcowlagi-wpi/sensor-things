%% Kalman Filter implementation for 1D navigation

clear variables; close all; clc

%%
% Load given data file; it must be in the same folder as this file.
% In practice, data are not available prerecorded like this example.
% Instead, data are received in real time.
% This data file loads the following variables with self-evident meanings:
% |time_stamps|, |accelerometer|, |gps_position|, and |gps_speed.|
load data_1D_navigation.mat

%% Error Covariances
%%
% Process noise is noise in accelerometer in this example,
% because we use acceleration data in the predictive model (kinematics).
% Process noise variance is a scalar here because 1D acceleration is
% measured.
% q	= (0.25)^2; %(m/s2)^2
% Q	= q;

q	= (0.05)^2; %(m/s2)^2
Q	= q;

%%
% Measurement noise is noise in the position and speed measurements.
% Therefore measurement noise covariance is a 2x2 square matrix with
% position and speed variance on diagonals.
R	= [1 0; 0 (0.5)^2]; % m^2 and (m/s)^2

%% Time Step and Interval of Interest
%%
% We have uniformly spaced time stamps, so calculating the difference
% between any two successive instants gives us $dt.$
dt	= time_stamps(2) - time_stamps(1);
%%
% Note the number of time stamps, which will be the number of iterations
% for which the Kalman filter will run.
n_t_pts = numel(time_stamps);

%% Predictive Model and Measurement Model
% Remember that the predictive model is 
%%
% $$ \underline{x}(k) = F\underline{x}(k-1) + G_1\underline{u}(k-1) + \underline{w}(k-1)$$
%%
% In this example $\underline{u}$ and $\underline{w}$ are scalars.
A	= [1 dt; 0 1];
B	= [0; dt];
G	= [0; dt];
%%
% The measurement model is $\underline{z} = C\underline{x}.$
C	= eye(2);

%% Initialization
%%
% Set initial values for state estimate $\underline{\hat{x}}$ and
% estimation error covariance $P$.
xHat	= [gps_position(1); gps_speed(1)];
P		= R;

[PSteadyState,~,~]	= idare(A', C', (G*Q*G'), R, [], []);
PSteadyStateStore	= reshape(PSteadyState, 4, 1);
LSteadyState		= PSteadyState * C' / (C*PSteadyState*C' + R);
%% 
% Memory allocation for recording the state estimate and error covariance
% over time. We will store the result of each iteration in the variables
% |xhat_store| and |P_store|. We also store the trace of the P matrix at
% each iteration in |P_trace_store|. Each column of these variables stores
% $\underline{\hat{x}}$ and $P$ at each iteration. 
storeXHat	= zeros(2, n_t_pts);
storeP		= zeros(4, n_t_pts);
storePTrace	= zeros(1, n_t_pts);
storeInnov	= zeros(2, n_t_pts);
%%
% Therefore, the first column stores the initial values. In general the
% $k^\mathrm{th}$ column of |xhat_store| stores $\underline{\hat{x}}(k-1)$
% and so on. The $(k-1)$ is because MATLAB's array numbering starts at 1,
% whereas our $k$ values (time stamps) start at 0.
storeXHat(:, 1)	= xHat;
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
for m1 = 1:(n_t_pts-1)
%%
% Prediction equations to get a preliminary estimate and error covariance.
	u	= accelerometer(m1);
	xHatMinus	= A*xHat + B*u;
	PMinus	= A*P*A' + G*Q*G';
%%	
% Compute Kalman gain. Note the use of |/| command instead of inverse.
	L	= (PMinus * C') / (C * PMinus * C' + R );
	
%% 
% Correction equations to get new estimate and error covariance at this
% iteration. Note the |eye(2)| command, which is a 2x2 identity matrix.
	z			= [gps_position(m1 + 1); gps_speed(m1 + 1)];
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

%% Innovation analysis and autocorrelation
innovCovar		= C * PSteadyState * C' + R;

nTimeStamps		= length(time_stamps);
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
fontSize_ = 20;

figure('Units','normalized', 'Position', [0.1 0.1 0.5*[1 1]]) 
subplot(211); hold on; plot(time_stamps, storeXHat(1,:), 'LineWidth', 2)
title('State estimate mean $\hat{x}_1$ (position)', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', fontSize_);
ylabel('$\hat{x}_1$ (m)', 'Interpreter', 'latex', 'FontSize', fontSize_);
ax = gca; ax.FontSize = fontSize_; ax.FontName = 'Times New Roman';

subplot(212); hold on; plot(time_stamps, storeXHat(2,:), 'LineWidth', 2)
title('State estimate mean $\hat{x}_2$ (speed)', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', fontSize_);
ylabel('$\hat{x}_2$ (m/s)', 'Interpreter', 'latex', 'FontSize', fontSize_);
ax = gca; ax.FontSize = fontSize_; ax.FontName = 'Times New Roman';

exportgraphics(gcf, 'kf_1Dnav_state_underestprocess.png', 'Resolution', 150)


figure;
plot(time_stamps, storePTrace, 'LineWidth', 2)
title('Trace of estimation error covariance $P$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', fontSize_);
ylabel('tr$(P)$, units N/A', 'Interpreter', 'latex', 'FontSize', fontSize_);
ax = gca; ax.FontSize = fontSize_; ax.FontName = 'Times New Roman';

figure;
plot(time_stamps, storeP, 'LineWidth', 2); hold on;
ax = gca;
ax.ColorOrderIndex = 1;
plot(time_stamps, kron(PSteadyStateStore, ones(1, length(time_stamps))), 'LineWidth', 2, 'LineStyle', '--')
ax = gca; ax.FontSize = fontSize_;
title('Elements of estimation error covariance $P$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', fontSize_);
ylabel('$p_{ij}$', 'Interpreter', 'latex', 'FontSize', fontSize_);
legend('$p_{11}$', '$p_{12}$', '$p_{21}$', '$p_{22}$', 'Interpreter', 'latex' )
ax = gca; ax.FontSize = fontSize_; ax.FontName = 'Times New Roman';


figure('Units','normalized', 'Position', [0.1 0.1 0.5*[1 1]])
plot(time_stamps, storeInnov, 'LineWidth', 2); hold on;
ax = gca; ax.ColorOrderIndex = 1;
plot(time_stamps, 2*sqrt(innovCovar(1,1))*ones(1, length(time_stamps)), ...
	'LineWidth', 2, 'LineStyle', '--'); ax.ColorOrderIndex = 1;
plot(time_stamps, -2*sqrt(innovCovar(1,1))*ones(1, length(time_stamps)), ...
	'LineWidth', 2, 'LineStyle', '--')
plot(time_stamps, 2*sqrt(innovCovar(2,2))*ones(1, length(time_stamps)), ...
	'LineWidth', 2, 'LineStyle', '--'); ax.ColorOrderIndex = 2;
plot(time_stamps, -2*sqrt(innovCovar(2,2))*ones(1, length(time_stamps)), ...
	'LineWidth', 2, 'LineStyle', '--')
% title('Innovations sequence', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', fontSize_);
ylabel('Innovation', 'Interpreter', 'latex', 'FontSize', fontSize_);
ax = gca; ax.FontSize = fontSize_; ax.FontName = 'Times New Roman';
exportgraphics(gcf, 'kf_1Dnav_innov_underestprocess.png', 'Resolution', 150)

figure('Units','normalized', 'Position', [0.1 0.2 0.5*[1 1]])
plot(time_stamps, innovAutoCorr, 'LineWidth', 2); hold on;
plot(time_stamps, twoSigmaBd, 'LineWidth', 2, 'LineStyle', '--' );
ax = gca; ax.FontSize = fontSize_;
ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot(time_stamps, -twoSigmaBd, 'LineWidth', 2, 'LineStyle', '--' )
% title('Innovations autocorrelation', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', fontSize_);
ylabel('Innovation', 'Interpreter', 'latex', 'FontSize', fontSize_);
ax = gca; ax.FontSize = fontSize_; ax.FontName = 'Times New Roman';
exportgraphics(gcf, 'kf_1Dnav_autocorr_underestprocess.png', 'Resolution', 150)

