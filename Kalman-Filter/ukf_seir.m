function ukf_seir()
% Unscented Kalman Filter implementation for SEIR contagion model

close all; clc

%% Problem setup and data
load data_seir.mat zMeas xTrue nTimeStamps timeStamps stdDevProcNoise stdDevMeasNoise

%----- Dimensions
nStates = 4;		% States are Susceptible, Infected, Exposed, Removed fractions
nMeas	= 1;		% Measurement is Infected fraction

nProcNoise	= 4;	% All state equations are noisy with equal variance and no correlation
nXAug		= nStates + nProcNoise + nMeas;

%----- SEIR model parameters
% Different diseases have different coefficients

SEIRmodelParameters.beta	= 0.5;
SEIRmodelParameters.alpha	= 0.1;
SEIRmodelParameters.delta	= 0.05;


%----- Time stamps
dt_	= timeStamps(2) - timeStamps(1);	% assuming uniformly spaced time stamps

C	= [0 0 1 0];

%% Error Covariances
QCont	= stdDevProcNoise^2*eye(4);					% covariance matrix of the continuous process noise
R		= stdDevMeasNoise^2 ;						% variance of the measurement noise

Q		= QCont*dt_;								% first-order approximation of discrete-time proc. noise covariance

%% Initialization
xHat		= [1 - zMeas(1); 0; zMeas(1); 0];		% this is model-specific; we assume Susceptible = 1 - Infected fraction; no one Exposed or Removed
P			= 1E-3*eye(4);								% initialize based on measurement covariance

storeXHat	= zeros(nStates, nTimeStamps);
storeP		= zeros(nStates^2, nTimeStamps);
storePTrace	= zeros(1, nTimeStamps);
storeInnov	= zeros(nMeas, nTimeStamps);
storeInnovCovar	= zeros(nMeas^2, nTimeStamps);

storeXHat(:, 1)		= xHat;
storePTrace(:, 1)	= trace(P);
storeP(:, 1)		= reshape(P, nStates^2, 1); 

%% UKF initialization and weights
xSigma				= zeros(nStates, 2*nXAug);
zSigma				= zeros(nMeas, 2*nXAug);
weightUT0			= 1/3;
weightUT			= (1 - weightUT0) / (2*nXAug);

%% Run UKF iterations
for m1 = 1:(nTimeStamps-1)

	%----- Generate sigma points and propagate through system equations
	PXAug		= blkdiag(P, Q, R);
	xAugSigma0	= [xHat; zeros(nProcNoise, 1); zeros(nMeas, 1)];			% sigma points of augmented state

	xSigma0		= rk4_step( timeStamps(m1), dt_, ...
		xAugSigma0(1:nStates), [], ...
		xAugSigma0(nStates+1 : nStates+nProcNoise), ...
		@seir_dynamics, SEIRmodelParameters );								% We use RK4 for one-step propagation, OK to use simpler Euler integration as well
	
	xAugSigma	= generate_sigma_points(xAugSigma0, PXAug);
	for m2 = 1:2*nXAug
		xSigma(:, m2)	= rk4_step( timeStamps(m1), dt_, ...
			xAugSigma(1:nStates, m2), [], ...								% State
			xAugSigma(nStates+1 : nStates+nProcNoise, m2), ...				% Process noise
			@seir_dynamics, SEIRmodelParameters );
	end
	
	%----- Predictive updates using sigma points weighted sum
	xHatMinus	= weightUT0*xSigma0 + weightUT*sum(xSigma, 2);
	PMinus		= weightUT0*(xSigma0 - xHatMinus)*(xSigma0 - xHatMinus)';
	for m2 = 1:2*nXAug		% this can be vectorized by reshaping the matrix into a column vector
		PMinus	= PMinus + ...
			weightUT*(xSigma(:, m2) - xHatMinus)*(xSigma(:, m2) - xHatMinus)';	% Note the transpose at the end
	end

	%----- Measurement model applied to propagated sigma points
	% We don't really need this for this particular example because the
	% measurement model is linear; including it for the sake of example
	zSigma0		= measurement_model(...
		xSigma0, xAugSigma0(nStates+nProcNoise+1 : nXAug));
	for m2 = 1:2*nXAug
		zSigma(:, m2)	= measurement_model(...
			xAugSigma(1:nStates, m2), ...									% State
			xAugSigma(nStates+nProcNoise+1 : nXAug, m2));					% Measurement noise
	end
	zHatMinus	= weightUT0*zSigma0 + weightUT*sum(zSigma, 2);

% 	zHatMinus	= C*xHatMinus;												% This is the usual when measurement model is linear

	%----- Error covariance and cross-covariance from sigma points
	PZZ			= weightUT0*(zSigma0 - zHatMinus)*(zSigma0 - zHatMinus)';		% We don't need to do this if meas. model is linear
	PXZ			= weightUT0*(xSigma0 - xHatMinus)*(zSigma0 - zHatMinus)';
	for m2 = 1:2*nXAug
		PZZ	= PZZ + ...
			weightUT*(zSigma(:, m2) - zHatMinus)*(zSigma(:, m2) - zHatMinus)';	
		PXZ	= PXZ + ...
			weightUT*(xSigma(:, m2) - xHatMinus)*(zSigma(:, m2) - zHatMinus)';	
	end

	%----- Compute Kalman gain. Note the use of |/| command instead of inverse.
% 	L			= (PMinus * C') / (C * PMinus * C' + R );					% This is the "usual" equation when meas. model is linear
	L			= PXZ / PZZ;
	
	%----- Measurement update
	z			= zMeas(m1 + 1);
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
			sum( shiftedInnov(1, m1+1:end).*storeInnov(1, m1+1:end) ) );
	twoSigmaBd(m1)			= 2/sqrt(nTimeStamps - m1); 
end

%% Figures of merit of the filter
format shortE

xTilde	= storeXHat - xTrue;		% this is the estimation error, which we can calculate because we know the ground truth in this example
sqEE	= xTilde(1, :).^2 + xTilde(2, :).^2 + xTilde(3, :).^2 + xTilde(4, :).^2;

fprintf('----- Mean of innovations autocorrelation (closer to zero is better): \n'); disp( mean(innovAutoCorr(2:end)) )
fprintf('----- Mean square estimation error (closer to zero is better): \n'); disp( mean(sqEE) )

disp(storePTrace(end))

%% Plot Results
figure;
xAxisTitles = {'Susceptible $x_1$', 'Exposed $x_2$', 'Infected $x_3$', 'Removed $x_4$'};
for m1 = 1:4
	subplot(2,2,m1);
	plot(timeStamps, storeXHat(m1,:), 'LineWidth', 2);
	hold on;
	plot(timeStamps, xTrue(m1, :), 'LineWidth',  2)
	make_nice_figures(gcf, gca, 18, [], 'Time', ...
		xAxisTitles{m1}, 'States', [0.1 0.24 0.5*[1 1]],[],[],[]);
	legend('Estimated', 'True')

end

% figure;
% plot(timeStamps, storePTrace, 'LineWidth', 2);
% make_nice_figures(gcf, gca, 18, [], 'Time', 'tr$(P)$', 'Trace', [0.15 0.35 0.5*[1 1]],[],[],[]);
% figure;
% plot(timeStamps, storeP, 'LineWidth', 2); hold on;
% make_nice_figures(gcf, gca, 18, [], 'Time', '$p_{ij}$', 'E.E. Covariance', [0.2 0.35 0.5*[1 1]],[],[],[]);

figure;
plot(timeStamps, storeInnov(1, :), 'LineWidth', 2); hold on;
ax = gca; ax.ColorOrderIndex = 1;
plot(timeStamps, 2*(storeInnovCovar(1,:).^0.5), ...
	'LineWidth', 2, 'LineStyle', '--'); ax.ColorOrderIndex = 1;
plot(timeStamps, -2*(storeInnovCovar(1,:).^0.5), ...
	'LineWidth', 2, 'LineStyle', '--')
make_nice_figures(gcf, gca, 18, 'Innovations Sequence', 'Time', ...
		'Innovation $\tilde{z}$', 'Innovations', [0.1 0.24 0.5*[1 1]],[],[],[]);

figure;
plot(timeStamps, innovAutoCorr, 'LineWidth', 2); hold on;
% plot(timeStamps, (2/sqrt(nTimeStamps))*ones(size(timeStamps)), 'LineWidth', 2, 'LineStyle', '--' );
% ax = gca; ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
% plot(timeStamps, (-2/sqrt(nTimeStamps))*ones(size(timeStamps)), 'LineWidth', 2, 'LineStyle', '--' )
make_nice_figures(gcf, gca, 18, 'Innovations Autocorrelation', 'Time', ...
		'Innovation a.c.', 'Innovation Autocorrelation', [0.1 0.24 0.5*[1 1]],[],[],[]);


%% Internal functions

	function zk_ = measurement_model(xk_, vk_)
		zk_ = [0 0 1 0]*xk_ + vk_;											% Linear model in this example
		% We don't really need the UT for the measurement update in this
		% specific example because the measurement model is linear
	end
	
	
	function xSigma_ = generate_sigma_points(xBar_, PX_)
		nX_		= numel(xBar_);												% In the UKF, the augmented state is different from the usual state
		S		= sqrt(nX_ / (1 - weightUT0)) * chol(PX_);
	
		xSigma_ = zeros(nX_, 2*nX_);
		for m10 = 1:nX_
			xSigma_(:, m10)			= xBar_ + S(:, m10);
			xSigma_(:, m10 + nX_)	= xBar_ - S(:, m10);
		end
	
	end


end