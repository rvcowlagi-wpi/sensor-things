function particle_propagation_ex1()

close all; clc;

aSys	= 1;
dt_		= 0.1;
dx_		= 0.5;

xHat	= 1;
P		= 0.1;
Q		= 0.1;
R		= 100;
G		= 1;

nParticles	= 1000;
particles	= xHat + sqrt(P)*randn(1, nParticles);
% weights_	= (1 / nParticles)*ones(1, nParticles);

% weights_*particles'

figure;
xPlot	= linspace(-2, 3, 1000);
% weight_ = (1 / nParticles);
for m1 = 1:6
	
	%----- These steps are for comparison / evaluation, not part
	% 	of the particle filter itself
	A	= exp(-aSys*dt_);
	xHat= A*xHat;
	P	= A*P*A' + G*Q*G';
	z	= measurement_model(xHat, sqrt(R)*randn);		% a measurement
	%-----------------------------------------------------------

	for m2 = 1:nParticles
		particles(m2)	= one_step_update(particles(m2), sqrt(Q)*randn);
% 		zHati			= measurement_model(particles(m2), 0);				% assuming zero mean measurement error
% 		weights_(m2)	= weights_(m2) * (1 / sqrt(2*pi*R)) * exp(-0.5*(z - zHati)^2 / R);
	end
% 	weights_	= weights_ / sum(weights_);
	
	particleDensity		= zeros(length(xPlot), 1);
	for m2 = 1:nParticles
		xIndex		= (xPlot >= particles(m2) - dx_ / 2) & (xPlot <= particles(m2) + dx_ / 2);
		particleDensity(xIndex) = particleDensity(xIndex) + (1 / nParticles)/dx_;
% 		particleDensity(xIndex) = particleDensity(xIndex) + weights_(m2)/dx_;
	end
	
	% Resampling needed
	
	yPlot	= normpdf(xPlot, xHat, sqrt(P));
	axString= ['$k = ' num2str(m1) '$'];
	subplot(2, 3, m1);
	plot(xPlot, particleDensity, '.'); hold on;
	plot(xPlot, yPlot, 'LineWidth', 3);
	make_nice_figures(gcf, gca, 14, axString, '$x$', '$p_X(x)$', [], [], [], [-2 3], [0 1.1])
end

	function xk_ = one_step_update(xkMinus1_, wkMinus1_)
		xk_	= exp(-aSys*dt_)*xkMinus1_ + wkMinus1_;
	end

	function zk_ = measurement_model(xk_, vk_)
		zk_	= xk_ + vk_;
	end

end
