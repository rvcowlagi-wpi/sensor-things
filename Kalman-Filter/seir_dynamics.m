function xDot = seir_dynamics(t_, x_, u_, v_, modelParameters)

	c1		= modelParameters.beta;
	c2		= modelParameters.alpha;
	c3		= modelParameters.delta;


	xDot	= zeros(4, 1);
	xDot(1)	= -c1 * x_(1) * x_(3);
	xDot(2)	= c1 * x_(1) * x_(3) - c2 * x_(2);
	xDot(3)	= c2 * x_(2) - c3 * x_(3);
	xDot(4)	= c3 * x_(3);

	xDot	= xDot + v_;

end