function sensor_latency()

close all; clc

%% Constants
a_	= 1;

%% Measurements
z1	= 0.8;
z2	= 1.15;
z3	= 0.4;

%% Marginal (constants)

c1	= (integral(@(x) (1 / a_).*x.*exp(-x.*z1), 0, a_))^(-1);
c2	= (integral(@(x) (c1 / a_).*(x.^2).*exp(-x.*(z1 + z2)), 0, a_))^(-1);
c3	= (integral(@(x) (c2*c1 / a_).*(x.^3).*exp(-x.*(z1 + z2 + z3)), 0, a_))^(-1);

%% LSEE
xLSEE1	= integral(@(x) x.*posterior_1(x, z1), 0, a_)
xLSEE2	= integral(@(x) x.*posterior_2(x, z2), 0, a_)
xLSEE3	= integral(@(x) x.*posterior_3(x, z3), 0, a_)

%% Linear LSEE with z1
% varZ		= 
% covXZ		= []
% xLinLSEE1	= 0.5*a_ + (covXZ / varZ)*(z1 - 1);	% E[Z] = 1

%% Sanity check; these should all be 1
fprintf('\n\n ---- Sanity Checks (everything should be 1) ---- \n')
integral(@prior_X, 0, a_)
integral(@(x) posterior_1(x, z1), 0, a_)
integral(@(x) posterior_2(x, z2), 0, a_)
integral(@(x) posterior_3(x, z3), 0, a_)


%% Plots
figure;
xPlot = linspace(0, a_, 100);

%----- Prior and posteriors
plot(xPlot, prior_X(xPlot), 'LineWidth', 2); hold on; ax = gca;
plot(xPlot, posterior_1(xPlot, z1), 'LineWidth', 2);
plot(xPlot, posterior_2(xPlot, z2), 'LineWidth', 2);
plot(xPlot, posterior_3(xPlot, z3), 'LineWidth', 2);

% 
% %----- Measurements
% zPlot	= linspace(0, 1.1*max(posterior_1(xPlot, z1)), 100);
% plot(z1*ones(size(xPlot)), zPlot, 'LineWidth', 2, 'LineStyle','--')
% plot(z2*ones(size(xPlot)), zPlot, 'LineWidth', 2, 'LineStyle','--')
% plot(z3*ones(size(xPlot)), zPlot, 'LineWidth', 2, 'LineStyle','--')
% 
% leg_ = legend('Prior', 'Posterior 1', 'Posterior 2', 'Posterior 3', '$z_1$', '$z_2$', '$z_3$');
% leg_.Interpreter = 'latex';

make_nice_figures(gcf, gca, 18, [], '$x$', '$f_X(x)$', 'PDF', [], [], [], []);



%% Functions

	function pdf_ = conditional_ZGivenX(x_, z_)
		pdf_ = x_.*exp( -x_.*z_);
	end

	
	function pdf_ = prior_X(x_)
		pdf_ = (1/a_)*ones(size(x_));
	end

	function pdf_ = posterior_1(x_, z_)
		pdf_ = c1*(1 / a_).*x_.*exp(-x_.*z_);
	end

	function pdf_ = posterior_2(x_, z_)
		pdf_ = c2*(c1 / a_).*(x_.^2).*exp(-x_.*(z1 + z_));
	end

	function pdf_ = posterior_3(x_, z_)
		pdf_ = c3*c2*(c1 / a_).*(x_.^3).*exp(-x_.*(z1 + z2 + z_));
	end
	

end