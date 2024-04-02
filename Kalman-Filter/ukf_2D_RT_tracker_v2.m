function ukf_2D_RT_tracker()

%% Unscented Kalman Filter implementation for 2D Range-Bearing Tracking
% Copy-pasted from 1D navigation; fix comments later

close all; clc

%%
load data_2D_RT_tracker.mat bearTrue rangeTrue timeStamps zBear zRange

nStates = 4;		% States are r, rDot, theta, thetaDot
nMeas	= 2;		% Measurements are r, theta
nProcNoise = 2;		% Say we consider noise in rDot and thetaDot processes
nXAug	= nStates + nProcNoise + nMeas;

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
% G	= [ [0; dt] zeros(2, 1); zeros(2, 1) [0; dt] ];
% Linearized G not needed in UKF
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

%%
storeXHat(:, 1)		= xHat;
storePTrace(:, 1)	= trace(P);
storeP(:, 1)		= reshape(P, nStates^2, 1); 

xSigma				= zeros(nStates, 2*nXAug);
zSigma				= zeros(nMeas, 2*nXAug);

const2UT			= 1E-2;		% alpha
const3UT			= 2;		% beta
const4UT			= 0;		% kappa
const1UT			= const2UT^2*(nXAug + const4UT)
lamUT0M				= const1UT / (nXAug + const1UT);
lamUT0C				= const1UT / (nXAug + const1UT) + (1 - const2UT^2 + const3UT);
lamUT				= 1 / ( 2*(nXAug + const1UT) );


%% Run Kalman Filter Iterations
for m1 = 1:(nTimeStamps-1)

	%----- Generate sigma points and propagate through system equations
	PXAug		= blkdiag(P, Q, R);
	xAugSigma0	= [xHat; zeros(nProcNoise, 1); zeros(nMeas, 1)];
	xSigma0		= one_step_update(...
		xAugSigma0(1:nStates), xAugSigma0(nStates+1 : nStates+nProcNoise));
	xAugSigma	= generate_sigma_points(xAugSigma0, PXAug);
	for m2 = 1:2*nXAug
		xSigma(:, m2)	= one_step_update(...
			xAugSigma(1:nStates, m2), ...									% State
			xAugSigma(nStates+1 : nStates+nProcNoise, m2) ...				% Process noise
			);
	end
	
	%----- Prediction equations to get a preliminary estimate and error covariance.
% 	u			= 0;
% 	A			= eye(nStates) + jacobianA(xHat, V)*dt;
% 	xHatMinus	= one_step_update(xkMinus1_);
% 	PMinus		= A*P*A' + G*Q*G';
	% The equations above are EKF update equations, for comparison

	%----- Predictive updates using sigma points weighted sum
	% Leave the process noise, take the states 
	%						-- Clemenza
	xHatMinus	= lamUT0M*xSigma0 + lamUT*sum(xSigma, 2);
	PMinus		= lamUT0C*(xSigma0 - xHatMinus)*(xSigma0 - xHatMinus)';
	for m2 = 1:2*nXAug
		PMinus	= PMinus + ...
			lamUT*(xSigma(:, m2) - xHatMinus)*(xSigma(:, m2) - xHatMinus)';	% Note the transpose at the end
	end

	%----- Measurement model applied to propagated sigma points
	zSigma0		= measurement_model(...
		xSigma0, xAugSigma0(nStates+nProcNoise+1 : nXAug));
	for m2 = 1:2*nXAug
		zSigma(:, m2)	= measurement_model(...
			xAugSigma(1:nStates, m2), ...									% State
			xAugSigma(nStates+nProcNoise+1 : nXAug, m2));					% Process noise
	end
	zHatMinus	= lamUT0M*zSigma0 + lamUT*sum(zSigma, 2);

	zHatMinus	= C*xHatMinus;

	%----- Error covariance and cross-covariance from sigma points
	PZZ			= lamUT0C*(zSigma0 - zHatMinus)*(zSigma0 - zHatMinus)';		% We don't need to do this if meas. model is linear
	PXZ			= lamUT0C*(xSigma0 - xHatMinus)*(zSigma0 - zHatMinus)';
	for m2 = 1:2*nXAug
		PZZ	= PZZ + ...
			lamUT*(zSigma(:, m2) - zHatMinus)*(zSigma(:, m2) - zHatMinus)';	
		PXZ	= PXZ + ...
			lamUT*(xSigma(:, m2) - xHatMinus)*(zSigma(:, m2) - zHatMinus)';	
	end

	%----- Compute Kalman gain. Note the use of |/| command instead of inverse.
	L			= (PMinus * C') / (C * PMinus * C' + R );					% This is the "usual" equation when meas. model is linear
% 	L			= PXZ / PZZ;
	
	%----- Measurement update
	z			= [zRange(m1 + 1); zBear(m1 + 1)];
	thisInnov	= z - zHatMinus;
	
	xHat		= xHatMinus + L*(thisInnov);
% 	P			= (eye(nStates) - L*C)*PMinus;								% This is the usual equation
	P			= PMinus - L*PZZ*L';										% Identical to the previous equation

	thisInnovCovar		= PZZ; 

	%----- Store results
	storeXHat(:, m1+1)		= xHat;
	storeP(:, m1+1)			= reshape(P, nStates^2, 1);
	storePTrace(:, m1+1)	= trace(P);
	storeInnov(:, m1+1)		= thisInnov;
	storeInnovCovar(:, m1+1)= reshape(thisInnovCovar, nMeas^2, 1);
end

%% Innovation autocorrelation
nTimeStamps		= length(timeStamps);

innovAutoCorr	= zeros(1, nTimeStamps);
shiftedInnov	= zeros(nMeas, nTimeStamps);
twoSigmaBd		= zeros(1, nTimeStamps);
for m1 = 1:nTimeStamps - 10
	shiftedInnov(:,m1+1:end)= storeInnov(:, 2:(nTimeStamps - m1 + 1));
	innovAutoCorr(m1)		= (1 / (nTimeStamps - m1) ) * (...
			sum( shiftedInnov(1, m1+1:end).*storeInnov(1, m1+1:end) ) + ...
			sum( shiftedInnov(2, m1+1:end).*storeInnov(2, m1+1:end) ) );
	twoSigmaBd(m1)			= 2/sqrt(nTimeStamps - m1); 
end

%% Figures of merit of the filter
format shortE

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
	'Range $\hat{x}_1 = r$ (km)', 'Range and Rate', [0.1 0.24 0.5*[1 1]],[],[],[]);

subplot(222); hold on; plot(timeStamps, storeXHat(2,:), 'LineWidth', 2)
make_nice_figures(gcf, gca, 18, [], 'Time (h)', ...
	'Range rate $\hat{x}_2 = \dot{r}$ (km/hr)', 'Range and Rate', [0.1 0.24 0.5*[1 1]],[],[],[]);

subplot(223); hold on; plot(timeStamps, storeXHat(3,:)*180/pi, 'LineWidth', 2);
hold on;
plot(timeStamps, bearTrue*180/pi, 'LineWidth',  2)
make_nice_figures(gcf, gca, 18, [], 'Time (h)', ...
	'Range $\hat{x}_3 = \theta$ (deg)', 'Bearing and Rate', [0.1 0.24 0.5*[1 1]],[],[],[]);

subplot(224); hold on; plot(timeStamps, storeXHat(4,:)*180/pi, 'LineWidth', 2)
make_nice_figures(gcf, gca, 18, [], 'Time (h)', ...
	'Range rate $\hat{x}_4 = \dot{\theta}$ (deg/hr)', 'Bearing and Rate', [0.1 0.24 0.5*[1 1]],[],[],[]);




fig2 = figure;
plot(timeStamps, storePTrace, 'LineWidth', 2);
make_nice_figures(gcf, gca, 18, [], 'Time (h)', 'tr$(P)$', 'Trace', [0.15 0.35 0.5*[1 1]],[],[],[]);

fig3 = figure;
plot(timeStamps, storeP, 'LineWidth', 2); hold on;
make_nice_figures(gcf, gca, 18, [], 'Time (h)', '$p_{ij}$', 'E.E. Covariance', [0.2 0.35 0.5*[1 1]],[],[],[]);

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
plot(timeStamps, -2*(storeInnovCovar(4, :).^0.5), ...
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
	function xk_ = one_step_update(xkMinus1_, wkMinus1_)
		xk_ = xkMinus1_ + ...
			system_dynamics_continuous([], xkMinus1_, wkMinus1_)*dt;		% First-order Euler approximation

		% Alternative to Euler approximation: RK4 (need to fix code to
		% accommodate w)
		% [~, xSim]	= ode45(@polar_kinematics_2D(t,x), [0 dt], xkMinus1_);
		% xk_	= xSim(end, :)';
	end

	function zk_ = measurement_model(xk_, vk_)
		zk_ = C*xk_ + vk_;													% Linear model in this example
		% We don't really need the UT for the measurement update in this
		% specific example because the measurement model is linear
	end
	
	function xDot = system_dynamics_continuous(t_, x_, w_)
		xDot	= zeros(4, 1);

		%----- Range-bearing kinematics
		% For a target moving at a constant velocity along
		% the 1st Cartesian axis
	
% 		xDot(1) = x_(2);
% 		xDot(3)	= x_(4);
		xDot(1)	= V*cos(x_(3));
		xDot(3)	= -V*sin(x_(3)) / x_(1);
		xDot(2)	= -V*x_(4)*sin(x_(3)) + w_(1);
		xDot(4)	= V*x_(2)*sin(x_(3)) / (x_(1)^2) - V*x_(4)*cos(x_(3)) / x_(1) + w_(2);
	end
	
	function xSigma_ = generate_sigma_points(xBar_, PX_)
		nX_		= numel(xBar_);												% In the UKF, the augmented state is different from the usual state
		S		= sqrt(const1UT + nX_) * chol(PX_);
	
		xSigma_ = zeros(nX_, 2*nX_);
		for m10 = 1:nX_
			xSigma_(:, m10)			= xBar_ + S(m10, :)';
			xSigma_(:, m10 + nX_)	= xBar_ - S(m10, :)';
		end
	
	end

end
