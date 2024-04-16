function static_range_ex2()

% Measurement error variance depends on range

%% Problem setup
nSensors= 1;
r0		= (0.1)^2;
r1		= 0.01;

%----- Prior
xHat	= [1.1; 1.25];
PXX		= [(0.1)^2 0; 0 (0.2)^2];


%% Various PDFs
	function pdf_ = likelihood_ZGivenX(x_, z_, sensorLocations)
		pdf_	= ones(1, size(z_, 2));
		for m01 = 1:nSensors
			range_	= norm(x_ - sensorLocations(:, m01));					% Range to target
			rMeas_	= meas_error_variance(range_);							% Measurement error variance as a function of range
			pdf_	= pdf_ .* ...											% We assume conditional independence among sensors
				(1 ./ (sqrt(2*pi*rMeas_) )) .* ...
				exp( -(z_ - range_).^2 / (2*rMeas_) );
		end
	end

	
	function pdf_ = prior_X(x_)
		for m01 = 1:size(x_, 2)
			pdf_(m01) = (1 / (2*pi*det(PXX))) * ...
				exp( -0.5*(x_(:, m01) - xHat)'*PXX*(x_(:, m01) - xHat));
		end
	end

	function pdf_ = joint_XY(x_, z_)
		pdf_ = prior_X(x_) .* likelihood_ZGivenX(x_, z_);
	end

	function pdf_ = marginal_Z(z_)
		pdf_ = zeros(size(z_));
		for k = 1:size(z_, 2)
			pdf_(k) = integral(@(x_) joint_XY(x_, z_(:, k)), 1, 1000);
		end
	end

	function r_	= meas_error_variance(range_)
		r_	= r0 + r1*range_;
	end

end