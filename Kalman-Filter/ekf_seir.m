function ekf_seir()
% Extended Kalman Filter implementation for SEIR contagion model

close all; clc

%% Problem setup and data
load data_seir.mat zMeas xTrue nTimeStamps timeStamps stdDevProcNoise stdDevMeasNoise

%----- Dimensions
nStates = 4;		% States are Susceptible, Infected, Exposed, Removed fractions
nMeas	= 1;		% Measurement is Infected fraction

nProcNoise	= 4;	% All state equations are noisy with equal variance and no correlation

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

%% Run EKF iterations
for m1 = 1:(nTimeStamps-1)

	%----- Predictive updates using system model and linearization (for covariance)
	xHatMinus	= rk4_step( timeStamps(m1), dt_, ...
		xHat, [], zeros(nProcNoise, 1), ...
		@seir_dynamics, SEIRmodelParameters );								% We use RK4 for updating the mean
	% Note that the EKF updates the mean by propagating zero noise (mean)
	% through the system dynamics; UKF and PF propagate noise as well

	A			= eye(nStates) + jacobianA_continuous(xHat)*dt_;
	PMinus		= A*P*A' + Q;												% G is identity, A is the Jacobian of the discrete-time model

	%----- Compute Kalman gain. Note the use of |/| command instead of inverse.
	PZZ			= (C * PMinus * C' + R );
	L			= (PMinus * C') / PZZ;										% This is the "usual" equation when meas. model is linear
	
	%----- Measurement update
	zHatMinus	= C*xHatMinus;												% This is the usual when measurement model is linear
	z			= zMeas(m1 + 1);
	thisInnov	= z - zHatMinus;
	
	xHat		= xHatMinus + L*(thisInnov);
	P			= (eye(nStates) - L*C)*PMinus;								% This is the usual equation
% 	P			= PMinus - L*PZZ*L';										% Identical to the previous equation

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
	function jacA = jacobianA_continuous(x_)
		jacA = zeros(4);

		jacA(1, :)	= -SEIRmodelParameters.beta*[x_(3) 0 x_(1) 0];
		jacA(2, :)	= SEIRmodelParameters.beta*[x_(3) 0 x_(1) 0] - ...
			SEIRmodelParameters.alpha*[0 1 0 0];
		jacA(3, :)	= [0 SEIRmodelParameters.alpha -SEIRmodelParameters.delta 0];
		jacA(4, :)	= [0 0 SEIRmodelParameters.delta 0];
	end


end