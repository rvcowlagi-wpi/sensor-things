function x_tPlusdt = rk4_step( t_, dt_, xt_, ut_, vt_, ...
	system_dynamics, modelParameters )

	% x, u, v, are state, control, and process noise terms

	a1	= 0.5;		a2	= 0.5;		a3	= 1;
	b1	= 0.5;		b2	= 0;		b3	= 0.5;
	b4	= 0;		b5	= 0;		b6	= 1;
	g1	= 1/6;		g2	= 1/3;		g3	= 1/3;		g4	= 1/6;

	k1	= dt_ * system_dynamics(t_, xt_, ut_, vt_, modelParameters);
	k2	= dt_ * system_dynamics(t_ + a1*dt_, xt_ + b1*k1, ut_, vt_, modelParameters);
	k3	= dt_ * system_dynamics(t_ + a2*dt_, xt_ + b2*k1 + b3*k2, ut_, vt_, modelParameters);
	k4	= dt_ * system_dynamics(t_ + a3*dt_, xt_ + b4*k1 + b5*k2 + b6*k3, ut_, vt_, modelParameters);

	x_tPlusdt = xt_ + g1*k1 + g2*k2 + g3*k3 + g4*k4;
end