function datagen_landmark_navigation()

%% Constants and parameters
close all; clc

linIncrement = 200;
t1	= linIncrement;
t2	= t1 + 10*pi;
t3	= t2 + linIncrement;
t4	= t3 + 10*pi;
t5	= t4 + linIncrement;
t6	= t5 + 10*pi;
t7	= t6 + linIncrement;

dt_			= 1E-2;					% s
timeStamps	= 0:dt_:t7;				% s
nTimeStamps	= length(timeStamps);

speedV		= 0.1;						% desired speed, m/s

stdDevProcNoise	= 0.03;

landmarks	= [10 10 -10 -10; 10 -10 10 -10]	% m
uControl	= zeros(1, nTimeStamps);
for m1 = 1:nTimeStamps
	if timeStamps(m1) <= t2 && timeStamps(m1) > t1
		uControl(m1)	= 0.1;
	elseif timeStamps(m1) > t3 && timeStamps(m1) <= t4
		uControl(m1)	= -0.1;
	elseif timeStamps(m1) > t5 && timeStamps(m1) <= t6
		uControl(m1)	= 0.1;
	end
end
u	= uControl + stdDevProcNoise*randn(1, nTimeStamps);

figure
plot(timeStamps, u)

x0	= [-12; -8; 0];
xSim= zeros(nTimeStamps, 3);
t	= 0;
xt	= x0;
xSim(1, :) = x0';

for m1 = 2:nTimeStamps
	xt	= rk4_step( t, xt, m1 - 1 );
	t	= timeStamps(m1);

	xSim(m1, :) = xt';
end

groundTruthPosition		= xSim(:, 1:2)';
groundTruthHeading		= xSim(:, 3)';

stdDevRange	= 0.1;			% m

nLandmarks	= size(landmarks, 2);
groundTruthRanges	= zeros( nLandmarks, nTimeStamps );
measuredRanges		= zeros( nLandmarks, nTimeStamps );
for m1 = 1:nLandmarks
	groundTruthRanges(m1, :)	= (...
		(xSim(:, 1) - landmarks(1, m1)).^2 + ...
		(xSim(:, 2) - landmarks(2, m1)).^2 ).^0.5;
	measuredRanges(m1, :)		= groundTruthRanges(m1, :) + ...
		stdDevRange*randn(1, nTimeStamps);
end

plot(timeStamps, groundTruthRanges(1, :)); hold on;
plot(timeStamps, measuredRanges(1,:), '.')

save assignment4_problem4.mat measuredRanges ...
	stdDevRange timeStamps nTimeStamps ...
	groundTruthRanges groundTruthPosition groundTruthHeading ...
	stdDevProcNoise landmarks u


figure;
plot( xSim(:, 1), xSim(:, 2), 'LineWidth', 2 );
hold on; axis equal
plot(landmarks(1, :), landmarks(2, :), 'x')
% plot( measuredRange.*cos(measuredBearing), measuredRange.*sin(measuredBearing), '.'); 
make_nice_figures(gcf, gca, 14, 'Orbit', '$y_1$ (km)', '$y_2$ (km)', [], [], [], [], [])

%%
	function xDot = planar_kinematics(t_, x_, tIndex_)
		ut_		= u(tIndex_);

		xDot(1, 1)	= speedV * cos(x_(3));
		xDot(2, 1)	= speedV * sin(x_(3));
		xDot(3, 1)	= ut_;
	end

	function x_tplus_dt = rk4_step( t_, xt_, tIndex_ )

		a1	= 0.5;		a2	= 0.5;		a3	= 1;
		b1	= 0.5;		b2	= 0;		b3	= 0.5;
		b4	= 0;		b5	= 0;		b6	= 1;
		g1	= 1/6;		g2	= 1/3;		g3	= 1/3;		g4	= 1/6;

		k1	= dt_ * planar_kinematics(t_, xt_, tIndex_);
		k2	= dt_ * planar_kinematics(t_ + a1*dt_, xt_ + b1*k1, tIndex_);
		k3	= dt_ * planar_kinematics(t_ + a2*dt_, xt_ + b2*k1 + b3*k2, tIndex_);
		k4	= dt_ * planar_kinematics(t_ + a3*dt_, xt_ + b4*k1 + b5*k2 + b6*k3, tIndex_);

		x_tplus_dt = xt_ + g1*k1 + g2*k2 + g3*k3 + g4*k4;
	end
end