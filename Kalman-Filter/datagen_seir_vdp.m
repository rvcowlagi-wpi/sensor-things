function datagen_seir_vdp()

close all; clc

%% SEIR

c1	= 0.5;
c2	= 0.1;
c3	= 0.05;
SEIRmodelParameters.beta	= c1;
SEIRmodelParameters.alpha	= c2;
SEIRmodelParameters.delta	= c3;

x0	= [0.99; 0; 0.01; 0];

Q	= (1E-2)^2;
R	= (1E-3)^2;

nTimeStamps = 1E3 + 1;
timeStamps	= linspace(0, 150, nTimeStamps);
dt_			= timeStamps(2) - timeStamps(1);

xSim		= zeros(nTimeStamps, 4);
xSim(1, :)	= x0';
x			= x0;
zSim		= zeros(nTimeStamps, 1);
zSim(1, :)	= x0(3) + sqrt(R)*randn;

for m1 = 1:(nTimeStamps - 1)

	x = rk4_step( timeStamps(m1), dt_, x, [], sqrt(Q)*randn(4, 1), ...
		@seir_dynamics, SEIRmodelParameters );
	xSim(m1 + 1, :)		= x';
	zSim(m1 + 1, :)		= x(3) + sqrt(R)*randn;

end

xTrue	= xSim';
zMeas	= zSim';

stdDevProcNoise = sqrt(Q);
stdDevMeasNoise = sqrt(R);

save data_seir.mat zMeas xTrue nTimeStamps timeStamps stdDevProcNoise stdDevMeasNoise



%%

	

	function xDot = vdp_dynamics(t_, x_, u_, v_)
		xDot	= zeros(2, 1);
		xDot(1) = x_(2);
		xDot(2) = mu_*(1 - x_1^2)*x_2 - x_1 + v_;
	end

end
