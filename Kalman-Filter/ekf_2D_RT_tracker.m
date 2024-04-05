%% Extended Kalman Filter implementation for 2D Range-Bearing Tracking
% Copy-pasted from 1D navigation; fix comments later

function ekf_2D_RT_tracker()

close all; clc

%%
load data_2D_RT_tracker.mat bearTrue rangeTrue timeStamps zBear zRange

nStates = 4;	% States are r, rDot, theta, thetaDot
nMeas	= 2;	% Measurements are r, theta

V		= 50;

%% Error Covariances
Q	= [(1E-2)^2 0; 0 (1E-2*pi/180)^2];
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
P		= diag([R(1,1), 0.1, R(2,2), 0.1]);
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
storeNormInnov	= zeros(1, nTimeStamps);

%%
storeXHat(:, 1)		= xHat;
storePTrace(:, 1)	= trace(P);
storeP(:, 1)		= reshape(P, nStates^2, 1); 

%% Run Kalman Filter Iterations
for m1 = 1:(nTimeStamps-1)

	%----- Prediction equations to get a preliminary estimate and error covariance.
	u			= 0;
	A			= eye(nStates) + jacobianA(xHat)*dt;

% 	xHatMinus	= xHat + polar_kinematics_2D([], xHat, V)*dt;				% First-order Euler approximation
% 	[~, xSim]	= ode45(@(t,x) polar_kinematics_2D(t,x, V), [0 dt], xHat);	% RK4
% 	xHatMinus	= xSim(end, :)';
	xHatMinus	= one_step_update(xHat, [0; 0]);
	PMinus		= A*P*A' + G*Q*G';

	%----- Compute Kalman gain. Note the use of |/| command instead of inverse.
	L			= (PMinus * C') / (C * PMinus * C' + R );
	
	%----- Measurement update
	z			= [zRange(m1 + 1); zBear(m1 + 1)];
	thisInnov	= z - C*xHatMinus;	

	thisInnovCovar		= C * PMinus * C' + R;
	thisNormalizedInnov = ( thisInnov' / thisInnovCovar ) * thisInnov;

	if m1 == 11
		disp(PMinus)
		disp(R)
		disp(z)
		disp(thisInnov)
		disp(thisInnovCovar)
		disp(thisNormalizedInnov)

		zC1 = [19; 0.6];
		zC2 = [18.82; 0.7];
		zC3	= [19.1; 0.56];

		zT1	= zC1 - C*xHatMinus;
		zT2	= zC2 - C*xHatMinus;
		zT3	= zC3 - C*xHatMinus;

		ni1 = ( zT1' / thisInnovCovar ) * zT1
		ni2 = ( zT2' / thisInnovCovar ) * zT2
		ni3 = ( zT3' / thisInnovCovar ) * zT3

		niSum	= thisNormalizedInnov + ni1 + ni3;
		bta0 = thisNormalizedInnov / niSum
		bta1 = ni1 / niSum
		bta3 = ni3 / niSum

		xHat0 = xHatMinus + L*(thisInnov);
		xHat1 = xHatMinus + L*(zT1);
		xHat3 = xHatMinus + L*(zT3);

		xHat0
		bta0*xHat0 + bta1*xHat1 + bta3*xHat3

	end

	xHat		= xHatMinus + L*(thisInnov);
	P			= (eye(nStates) - L*C)*PMinus;
	
	%----- Store results
	storeXHat(:, m1+1)		= xHat;
	storeP(:, m1+1)			= reshape(P, nStates^2, 1);
	storePTrace(:, m1+1)	= trace(P);
	storeInnov(:, m1+1)		= thisInnov;
	storeInnovCovar(:, m1+1)= reshape(thisInnovCovar, nMeas^2, 1);
	storeNormInnov(:, m1+1)	= thisNormalizedInnov;
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

%% Figures of merit of the filter
% format shortE

mseeRange	= mean( (storeXHat(1,:) - rangeTrue).^2 );
mseeBrng	= mean( (storeXHat(3,:) - bearTrue).^2 );

fprintf('----- Mean of innovations autocorrelation (closer to zero is better): \n'); disp( mean(innovAutoCorr(2:end)) )
fprintf('----- Mean square estimation error in range (closer to zero is better): \n'); disp( mseeRange )
fprintf('----- Mean square estimation error in bearing (closer to zero is better): \n'); disp( mseeBrng )


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


figure;
histogram(storeNormInnov, 'Normalization','pdf', 'NumBins', 10); hold on;
xPlot = linspace(0, 15, 100);
plot( xPlot, chi2pdf(xPlot, nMeas), 'LineWidth', 3);
make_nice_figures(gcf, gca, 14, [], 'Normalized innovation', 'PDF', 'Normalized Innovation', [],[],[],[]);

figure;
plot(xPlot, chi2cdf(xPlot, nMeas), 'LineWidth', 3); grid on;
make_nice_figures(gcf, gca, 14, [], 'Normalized innovation', 'CDF', 'Chi-square CDF', [],[],[],[]);


%%
	function xDot = system_dynamics_continuous(t_, x_, w_)
		xDot	= zeros(4, 1);
	
	% 	xDot(1) = x_(2);
	% 	xDot(3)	= x_(4);
		xDot(1)	= V*cos(x_(3));
		xDot(3)	= -V*sin(x_(3)) / x_(1);
		xDot(2)	= -V*x_(4)*sin(x_(3)) + w_(1);
		xDot(4)	= V*x_(2)*sin(x_(3)) / (x_(1)^2) - V*x_(4)*cos(x_(3)) / x_(1) + w_(2);
	end
	
	function jacA = jacobianA(x_)
		jacA = zeros(4);
		jacA(1, :)	= [0 1 0 0];
		jacA(3, :)	= [0 0 0 1];
		jacA(2, :)	= V*[0 0 -x_(4)*cos(x_(3)) -sin(x_(3))];
		jacA(4, :)	= V*[ ...
			-2*x_(2)*sin(x_(3))/(x_(1)^3) + x_(4)*cos(x_(3))/(x_(1)^2) ...
			sin(x_(3))/(x_(1)^2) ...
			x_(2)*cos(x_(3))/(x_(1)^2) + x_(4)*sin(x_(3))/x_(1) ...
			-cos(x_(3)) / x_(1)];
	end

	function xk_ = one_step_update(xkMinus1_, wkMinus1_)
		xk_ = xkMinus1_ + ...
			system_dynamics_continuous([], xkMinus1_, wkMinus1_)*dt;		% First-order Euler approximation

		% Alternative to Euler approximation: RK4 (need to fix code to
		% accommodate w)
		% [~, xSim]	= ode45(@polar_kinematics_2D(t,x), [0 dt], xkMinus1_);
		% xk_	= xSim(end, :)';
	end

end