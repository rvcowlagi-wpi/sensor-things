function x_ = asg5_target_step(t, xt_)

	x_ = asg5_rk4_step( t_, xt_ );

	function xDot = ag5_planar_kinematics(t_, x_)
		ut_		= [];

		speedV	= 0.1;

		xDot(1, 1)	= speedV * cos(x_(3));
		xDot(2, 1)	= speedV * sin(x_(3));
		xDot(3, 1)	= ut_;
	end

	function x_ = asg5_rk4_step( t_, xt_ )

		a1	= 0.5;		a2	= 0.5;		a3	= 1;
		b1	= 0.5;		b2	= 0;		b3	= 0.5;
		b4	= 0;		b5	= 0;		b6	= 1;
		g1	= 1/6;		g2	= 1/3;		g3	= 1/3;		g4	= 1/6;

		k1	= dt_ * ag5_planar_kinematics(t_, xt_, tIndex_);
		k2	= dt_ * ag5_planar_kinematics(t_ + a1*dt_, xt_ + b1*k1, tIndex_);
		k3	= dt_ * ag5_planar_kinematics(t_ + a2*dt_, xt_ + b2*k1 + b3*k2, tIndex_);
		k4	= dt_ * ag5_planar_kinematics(t_ + a3*dt_, xt_ + b4*k1 + b5*k2 + b6*k3, tIndex_);

		x_	= xt_ + g1*k1 + g2*k2 + g3*k3 + g4*k4;
	end


end
