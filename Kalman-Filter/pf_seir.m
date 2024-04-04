function pf_seir()
% Particle Filter implementation for SEIR contagion model

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

%% Error Covariances
QCont	= stdDevProcNoise^2*eye(4);					% covariance matrix of the continuous process noise
R		= (1E-2)^2;						% variance of the measurement noise

Q		= 1E-4*eye(4);								% arbitrary tuning to solve the "sample impoverishment" problem
sqrtQ	= chol(Q);
%% Initialization
xHat		= [1 - zMeas(1); 0; zMeas(1); 0];		% this is model-specific; we assume Susceptible = 1 - Infected fraction; no one Exposed or Removed
P			= 1E-3*eye(4);							% initialize based on measurement covariance

storeXHat	= zeros(nStates, nTimeStamps);
storeP		= zeros(nStates^2, nTimeStamps);
storePTrace	= zeros(1, nTimeStamps);

storeXHat(:, 1)		= xHat;
storePTrace(:, 1)	= trace(P);
storeP(:, 1)		= reshape(P, nStates^2, 1);

%% Initial particles

nParticles	= 1000;
xParticles	= xHat + chol(P)*randn(nStates, nParticles);							% Sample from the prior (assumed Gaussian)
particleWeights = (1 / nParticles) * ones(1, nParticles);					% Initialize with uniform weights
xResampled  = zeros(nStates, nParticles);
weightScaling = 100;

%% Run EKF iterations
for m1 = 1:(nTimeStamps-1)
	procNoiseSamples	= sqrtQ * randn(nProcNoise, nParticles);
	measNoiseSamples	= sqrt(R) * randn(nMeas, nParticles);

	%----- New measurement
	z	= zMeas(m1 + 1);

	%----- For each particle...
	for m2 = 1:nParticles
		xNext	= rk4_step( timeStamps(m1), dt_, ...
			xParticles(:, m2), [], procNoiseSamples(:, m2), ...
			@seir_dynamics, SEIRmodelParameters );							% propagate through system dynamics
		if any( isnan(xNext) )
			xNext
			xParticles(:, m2)
			procNoiseSamples(:, m2)
			fprintf('Huh?\n')
		end
		xParticles(:, m2) = xNext;

		zHati	= measurement_model(xParticles(:, 2), measNoiseSamples(:, m2));										% then through the measurement model
		particleWeights(m2)	= particleWeights(m2) * weightScaling * ...
			(1 / sqrt(2*pi*R)) * exp(-0.5*(z - zHati)^2 / R);				% Update weights proportional to likelihood
	end
% 	sum(particleWeights)
	particleWeights = particleWeights / sum(particleWeights);				% Normalize the weights to add up to 1
	
	

	%----- Resample particles proportional to weights and reset weights
	% Resampling means just picking particles from the original set with
	% replacement, with probability proportional to its weight

	fprintf('Iteration number %i, max. particle weight is %f: \n', m1, max(particleWeights));
	cumulWeights	= cumsum(particleWeights);
	for m2 = 1:nParticles
		a_		= rand;										% uniform random number in [0, 1]
		pIndex	= find( (a_ <= cumulWeights), 1, "first");	% particle number for which cumulative weight first exceeds the previous random number
		xResampled(:, m2)	= xParticles(:, pIndex);
	end
	xParticles = xResampled;



	%----- Store results
	% Note: the PF doesn't care about mean and covariance. It tries to
	% approximate the entire posterior distribution of the state. But we
	% generally keep track of the mean and covarianace summary statistics,
	% which can be easily computed from the particles. Among other things,
	% we can use it for plotting results.

	%----- Predictive updates using system model and linearization (for covariance)
	xHat	= (1 / nParticles) * sum(xParticles, 2);
	P		= zeros(nStates);
	for m2 = 1:nParticles
		P	= P + (xParticles(:, m2) - xHat) * ((xParticles(:, m2) - xHat)'); 
	end
	P		= (1 / (nParticles - 1) ) * P;
	
	%----- Store results
	storeXHat(:, m1+1)		= xHat;
	storeP(:, m1+1)			= reshape(P, nStates^2, 1);
	storePTrace(:, m1+1)	= trace(P);
end

%% Figures of merit of the filter
% format shortE

xTilde	= storeXHat - xTrue;		% this is the estimation error, which we can calculate because we know the ground truth in this example
sqEE	= xTilde(1, :).^2 + xTilde(2, :).^2 + xTilde(3, :).^2 + xTilde(4, :).^2;

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

figure;
plot(timeStamps, storePTrace, 'LineWidth', 2);
make_nice_figures(gcf, gca, 18, [], 'Time', 'tr$(P)$', 'Trace', [0.15 0.35 0.5*[1 1]],[],[],[]);
figure;
plot(timeStamps, storeP, 'LineWidth', 2); hold on;
make_nice_figures(gcf, gca, 18, [], 'Time', '$p_{ij}$', 'E.E. Covariance', [0.2 0.35 0.5*[1 1]],[],[],[]);


%% Internal functions	
	function zk_ = measurement_model(xk_, vk_)
		zk_ = [0 0 1 0]*xk_ + vk_;											% Linear model in this example
		% We don't really need the UT for the measurement update in this
		% specific example because the measurement model is linear
	end


end