function BREx_speed_radar_2()

close all; clc

% Measurement error has speed-dependent variance

r_	= 5^2; % (km/hr)^2
xSS	= linspace(0, 400, 100);
ySS	= linspace(-50, 400, 100);

infty = 1000;	% Practical "infinity"

[xMesh, yMesh]	= meshgrid(xSS, ySS);

fXPrior			= prior_X(xSS);
fYGivenX		= conditional_YGivenX(xSS(100), ySS);
fJoint			= joint_XY(xMesh, yMesh);

zStar			= 50;
xStar			= 60;
integral(@(x) posterior_XGivenY(x, zStar), xStar, infty)
% This is P[ X >= xBar | y = yMeas ]


% Sanity check; these should all be 1
fprintf('\n\n ---- Sanity Check ---- \n')
integral(@prior_X, 0, infty)
integral(@(y) conditional_YGivenX(xSS(50), y), 0, infty)
integral(@marginal_Y, -infty, infty)





fontSize_ = 40;
figure('units', 'normalized', 'OuterPosition', [0.05 0.05 0.8 0.9]);
surf(xSS, ySS, fJoint, 'EdgeColor','none')
xlabel('$x$', 'FontName', 'Times New Roman', ...
	'FontSize', fontSize_, 'FontWeight', 'bold', 'interpreter', 'latex');
ylabel('$y$', 'FontName', 'Times New Roman', ...
	'FontSize', fontSize_, 'FontWeight', 'bold', 'interpreter', 'latex');
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = fontSize_;


	function pdf_ = conditional_YGivenX(x_, y_)
		pdf_ = (1 ./ (sqrt(2*pi) .* (x_/10) ) ) .* ...
			exp( -(y_ - x_).^2 ./ (2 .* (x_.^2 / 100) ) );
	end

	
	function pdf_ = prior_X(x_)
		pdf_ = (1/50) * exp(-x_ / 50);
	end

	function pdf_ = joint_XY(x_, y_)
		pdf_ = prior_X(x_) .* conditional_YGivenX(x_, y_);
	end

	function pdf_ = marginal_Y(y_)
		pdf_ = zeros(size(y_));
		for k = 1:numel(y_)
			pdf_(k) = integral(@(x_) joint_XY(x_, y_(k)), 0, infty);
		end
	end

	function pdf_ = posterior_XGivenY(x_, y_)
		pdf_ = zeros(size(x_));
		for k = 1:numel(x_)
			pdf_(k) = joint_XY(x_(k), y_) ./ marginal_Y(y_);
		end
	end

end