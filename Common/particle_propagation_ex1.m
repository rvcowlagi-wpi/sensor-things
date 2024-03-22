function particle_propagation_ex1()

close all; clc;

aSys	= 1;
dt_		= 0.1;
dx_		= 0.01;

xHat	= 1;
P		= 0.1;
Q		= 0.1;
G		= 1;

nParticles	= 1000;
particles	= xHat + sqrt(P)*randn(1, nParticles);

figure;
xPlot	= linspace(-2, 3, 1000);
weight_ = (1 / nParticles);
for m1 = 1:6

	A	= exp(-aSys*dt_);

	xHat= A*xHat;
	P	= A*P*A' + G*Q*G';

	for m2 = 1:nParticles
		particlePropMean= one_step_update(particles(m2), 0);
		particles(m2)	= particlePropMean + sqrt(Q)*randn;
	end
	
	particleDensity		= zeros(length(xPlot), 1);
	for m2 = 1:nParticles
		xIndex		= (xPlot >= particles(m2)) & (xPlot <= particles(m2) + dx_);
		particleDensity = particleDensity + weight_* xPlot(xIndex);
	end
	
	
	yPlot	= normpdf(xPlot, xHat, sqrt(P));
	axString= ['$k = ' num2str(m1) '$'];
	subplot(2, 3, m1); plot(xPlot, yPlot, 'LineWidth', 3); hold on;
	plot(xPlot, particleDensity, 'LineWidth', 2)
	make_nice_figures(gcf, gca, 14, axString, '$x$', '$p_X(x)$', [], [], [], [-2 3], [0 1.1])
end

	function xk_ = one_step_update(xkMinus1_, wkMinus1_)
		xk_	= exp(-aSys*dt_)*xkMinus1_ + wkMinus1_;
	end

end