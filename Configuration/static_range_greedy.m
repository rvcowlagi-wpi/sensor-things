clear variables; close all; clc


%% Ground truth
xTrue	= [1; 1];

%% Prior
xHat	= [1.1; 1.25];
PXX		= [(0.1)^2 0; 0 (0.2)^2];

%% Mutual information

sensorMinRange	= 0.3;
sensorMaxRange	= 2;
rMeas			= (0.1)^2;

nPoints = 100;
x1F = (linspace(-1, 3.5, nPoints))';
x2F = linspace(-1, 3.5, nPoints);


nSensors		= 5;
sensorLocations = zeros(2, nSensors);
nSensorsPlaced	= 0;


figure;
plot(xHat(1), xHat(2), 'x', 'MarkerSize', 10); hold on; axis equal
xPlot	= linspace( (xHat(1) - (1-1E-6)*sqrt(PXX(1,1))), (xHat(1) + (1-1E-6)*sqrt(PXX(1,1))), 100 );
yPlot1	= sqrt(PXX(2,2)) * (1 - (xPlot - xHat(1)).^2 / PXX(1,1) ).^0.5 + xHat(2);
yPlot2	= -sqrt(PXX(2,2)) * (1 - (xPlot - xHat(1)).^2 / PXX(1,1) ).^0.5 + xHat(2);

plot(xPlot, yPlot1, 'LineWidth', 2);
ax = gca; ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot(xPlot, yPlot2, 'LineWidth', 2);

make_nice_figures(gcf, gca, 14, [], '$x_1$', '$x_2$', 'Mutual Information', [], [], [-1 3.5], [-1 3.5])

while (nSensorsPlaced < nSensors)
	maxMI = 0;
	maxMILocation = [];
	
	MI = zeros(nPoints);
	for m1 = 1:nPoints
		for m2 = 1:nPoints

			C = zeros(nSensorsPlaced + 1, 2);

			%----- These are the sensors already placed
			for m3 = 1:nSensorsPlaced
				thisSensorLocation = sensorLocations(:, m3);
				expRange = norm(thisSensorLocation - xHat);
				C(m3 + 1, :) = -(1 / expRange) * ...
					[(xHat(1) - thisSensorLocation(1)), (xHat(2) - thisSensorLocation(2))];
			end
		
			%----- This is the sensor to be placed in this iteration
			expRange	= norm([x1F(m1); x2F(m2)] - xHat);
			C(1, :)		= (1 / expRange)* ...
				[(xHat(1) - x1F(m1)), (xHat(2) - x2F(m2)) ];

			PZZ	= C * PXX * C' + rMeas*eye(nSensorsPlaced + 1);
			PXZ	= PXX*C';
			PJoint = [PXX PXZ; PXZ' PZZ];
	
			if (expRange <= sensorMaxRange) && (expRange >= sensorMinRange)
				MI(m2, m1) = 0.5*log( det(PXX) * det(PZZ) / det(PJoint) );
			end
	
			if MI(m2, m1) >= maxMI
				maxMI			= MI(m2, m1);
				maxMILocation	= [x1F(m1); x2F(m2)];
			end
		end
	end

	nSensorsPlaced = nSensorsPlaced + 1;
	sensorLocations(:, nSensorsPlaced) = maxMILocation;

	zMin	= min(MI(:));
	zMax	= max(MI(:));
	
	hdl1 = surf(ax, x1F, x2F, MI, 'LineStyle','none');
	if nSensorsPlaced > 1
		clim([1 zMax]);
	end
	colorbar; view(2); 
	hdl2 = plot3(maxMILocation(1), maxMILocation(2), zMax, '.', 'MarkerSize', 20);

	for m3 = 1:nSensorsPlaced - 1
		hdl3 = plot3(sensorLocations(1, m3), sensorLocations(2, m3), zMax, 'o', 'MarkerSize', 10);
	end

	filename_ = ['mi_example1_n' num2str(nSensorsPlaced) '.png'];
	exportgraphics(gca, filename_, 'Resolution', 150)

	delete(hdl1)
	delete(hdl2)
	if nSensorsPlaced > 1
		delete(hdl3)
	end

end







% 
