clear variables; close all; clc

priorXHat1	= [3; 5];
priorXHat2	= [7; 5];

priorP1		= [9 5; 5 9];
priorP2		= [25 -10; -10 9];

sensorLocations = [0 0; 5 0; 10 0]';
nSensors	= size(sensorLocations, 2);

R1			= [100 0; 0 (10*pi/180)^2];
R2			= [25 0; 0 (5*pi/180)^2];
R3			= [100 0; 0 (10*pi/180)^2];


figure;
for m1 = 1:nSensors
	plot(sensorLocations(1, m1), sensorLocations(2, m1), '.', 'MarkerSize', 20);
	hold on;
	ax = gca; ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
end
plot(priorXHat1(1), priorXHat1(2), 'x', 'MarkerSize', 10)
plot(priorXHat2(1), priorXHat2(2), 'x', 'MarkerSize', 10)

make_nice_figures(gcf, gca, 18, [], '$y_1$', '$y_2$', ...
	'Sensor FOV Configuration', [], [], [-1 11], [-1 11] )


S = (inv(priorP1));

x1Plot = linspace(-1, 11, 1E3);
x2Plot = linspace(-1, 11, 1E3);
for m1 = 1:1E3
	for m2 = 1:1E3
		x_ = [x1Plot(m1); x2Plot(m2)];
		if abs((x_ - priorXHat1)' * S * (x_ - priorXHat1) - 1) <= 5E-3
			plot(x_(1), x_(2), '.');
			ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
		end
	end
end

S = (inv(priorP2));
ax.ColorOrderIndex = ax.ColorOrderIndex + 1;
for m1 = 1:1E3
	for m2 = 1:1E3
		x_ = [x1Plot(m1); x2Plot(m2)];
		if abs((x_ - priorXHat2)' * S * (x_ - priorXHat2) - 1) <= 5E-3
			plot(x_(1), x_(2), '.');
			ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
		end
	end
end

%% Mutual information
%----- Covariances calculated approximately by linearization

Configuration: