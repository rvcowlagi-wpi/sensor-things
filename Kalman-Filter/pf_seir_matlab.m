function pf_seir_matlab()

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
R		= (1E-3)^2;						% variance of the measurement noise
Q		= (1E-4)*eye(4);				% arbitrary tuning to solve the "sample impoverishment" problem
sqrtQ	= chol(Q);

%% MATLAB Particle Filter
xHat		= [1 - zMeas(1); 0; zMeas(1); 0];

seirPF		= particleFilter(@seir_state_update_internal, @seir_likelihood_internal);
initialize(seirPF, 1000, xHat, eye(nStates));
seirPF.StateEstimationMethod	= 'mean';
seirPF.ResamplingMethod			= 'systematic';

storeXHat	= zeros(nStates, nTimeStamps);
for m1 = 1:nTimeStamps
    storeXHat(:, m1+1)		= correct(seirPF, zMeas(m1 + 1));
    predict(seirPF);
end

%% Internal functions
	function particles_ = seir_state_update_internal(particles_)

		[~, nParticles]		= size(particles_);
		procNoiseSamples	= sqrtQ * randn(nProcNoise, nParticles);

		for m2 = 1:nParticles
			xNext	= rk4_step( timeStamps(m1), dt_, ...
				particles_(:, m2), [], procNoiseSamples(:, m2), ...
				@seir_dynamics, SEIRmodelParameters );						% propagate through system dynamics
			particles_(:, m2) = xNext;
		end

	end

	function likelihood_ = seir_likelihood_internal(particles_, z_)

		[~, nParticles]		= size(particles_);
		likelihood_			= zeros(1, nParticles);

		for m2 = 1:nParticles
			zHati			= measurement_model(particles_(:, m2), 0);
			likelihood_(m2) = (1 / sqrt(2*pi*R)) * exp(-0.5*(z_ - zHati)^2 / R);
		end
	end

	function zk_ = measurement_model(xk_, vk_)
		zk_ = [0 0 1 0]*xk_ + vk_;											% Linear model in this example
		% We don't really need the UT for the measurement update in this
		% specific example because the measurement model is linear
	end

end
