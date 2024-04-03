function datagen_seir_vdp()

close all; clc

%%


%%

	function xDot = seir_dynamics(t_, x_, v_)
		xDot	= zeros(4, 1);
		xDot(1)	= ;

	end

	function xDot = vdp_dynamics(t_, x_, u_, v_)
		xDot	= zeros(2, 1);
		xDot(1) = x_(2);
		xDot(2) = mu_*(1 - x_1^2)*x_2 - x_1 + v_;
	end

end
