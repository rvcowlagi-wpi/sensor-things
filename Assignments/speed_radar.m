function speed_radar()

close all; clc

%% Constants
lambda_ = 1/50;			% (km/hr)^-1, inv. mean of prior
alpha_	= 0.1;			% std. dev. is alpha_*x
infty_	= 1000;			% Practically infinite

%% Posterior

z1		= 40;			% km/hr
x1		= 50;			% km/hr
ProbX	= integral(@(x) posterior_X(x, z1), x1, infty_);					% P[X >= 50]

extAP	= roots([lambda_ 0 100*z1 -100*z1^2]);
xMAP	= real(extAP(3))
xLSEE	= integral(@(x) x.*posterior_X(x, z1), 0, infty_)
xLinLSEE= (1 / lambda_) + ( (1 / lambda_^2) / (alpha_/lambda_ + 1 / lambda_^2) ) * (z1 - (1 / lambda_))


%% Plots
figure;
xPlot	= linspace(0, 100, 1000);
plot(xPlot, prior_X(xPlot), 'LineWidth', 2); hold on;
plot(xPlot, posterior_X(xPlot, z1), 'LineWidth', 2);
zPlot	= linspace(0, 1.1*max(posterior_X(xPlot, z1)), 1000);
plot(z1*ones(size(xPlot)), zPlot, 'LineWidth', 2, 'LineStyle','--')
plot(xMAP*ones(size(xPlot)), zPlot, 'LineWidth', 2, 'LineStyle','--')
plot(xLSEE*ones(size(xPlot)), zPlot, 'LineWidth', 2, 'LineStyle','--')
plot(xLinLSEE*ones(size(xPlot)), zPlot, 'LineWidth', 2, 'LineStyle','--')
make_nice_figures(gcf, gca, 18, [], '$x$', '$f_X(x)$', 'PDF', [], [], [], []);
legend('Prior', 'Posterior', 'Measurement', 'MAP', 'LSEE', 'Linear LSEE')



% % Measurement error has speed-dependent variance
% 
% xSS	= linspace(0, 400, 100);
% zSS	= linspace(-50, 400, 100);
% 
% infty = 1000;	% Practical "infinity"
% 
% [xMesh, yMesh]	= meshgrid(xSS, zSS);
% 
% fXPrior			= prior_X(xSS);
% fYGivenX		= conditional_YGivenX(xSS(100), zSS);
% fJoint			= joint_XY(xMesh, yMesh);
% 
% zStar			= 50;
% xStar			= 60;
% integral(@(x) posterior_XGivenY(x, zStar), xStar, infty)
% % This is P[ X >= xBar | y = yMeas ]


%% Sanity check; these should all be 1
fprintf('\n\n ---- Sanity Check ---- \n')
integral(@prior_X, 0, infty_)
integral(@(y) conditional_ZGivenX(x1, y), 0, infty_)
integral(@marginal_Y, -infty_, infty_)



	function pdf_ = conditional_ZGivenX(x_, z_)
		pdf_ = (1 ./ (sqrt(2*pi) .* (alpha_*x_)) ) .* ...
			exp( -0.5* ( (z_ - x_) ./ (alpha_*x_) ).^2 ); 
	end

	
	function pdf_ = prior_X(x_)
		pdf_ = lambda_ * exp(-lambda_*x_);
	end

	function pdf_ = joint_XY(x_, z_)
		pdf_ = prior_X(x_) .* conditional_ZGivenX(x_, z_);
	end

	function pdf_ = marginal_Y(z_)
		pdf_ = zeros(size(z_));
		for k = 1:numel(z_)
			pdf_(k) = integral(@(x_) joint_XY(x_, z_(k)), 0, infty_);
		end
	end

	function pdf_ = posterior_X(x_, z_)
		pdf_ = zeros(size(x_));
		for k = 1:numel(x_)
			pdf_(k) = joint_XY(x_(k), z_) ./ marginal_Y(z_);
		end
	end

end