clear variables; close all; clc


%% Ground truth
xTrue	= [1; 1];

%% Prior
xHat	= [1.1; 1.25];
PXX		= [(0.1)^2 0; 0 (0.2)^2];

%% Plot
figure;
plot(xHat(1), xHat(2), 'x', 'MarkerSize', 10); hold on; axis equal
xPlot	= linspace( (xHat(1) - (1-1E-6)*sqrt(PXX(1,1))), (xHat(1) + (1-1E-6)*sqrt(PXX(1,1))), 100 );
yPlot1	= sqrt(PXX(2,2)) * (1 - (xPlot - xHat(1)).^2 / PXX(1,1) ).^0.5 + xHat(2);
yPlot2	= -sqrt(PXX(2,2)) * (1 - (xPlot - xHat(1)).^2 / PXX(1,1) ).^0.5 + xHat(2);

plot(xPlot, yPlot1, 'LineWidth', 2);
ax = gca; ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot(xPlot, yPlot2, 'LineWidth', 2);

make_nice_figures(gcf, gca, 14, [], '$x_1$', '$x_2$', 'Fisher Information', [], [], [-1 3.5], [-1 3.5])


%% Calculate FIM contours for one sensor

sensorMinRange	= 0.3;
sensorMaxRange	= 2;
rMeas			= (0.1)^2;

% sensor2			= [2.9; 1.0];
% sensor2			= [2.25; 2.25];
sensor2			= [1.6; 1.6];

nPoints = 100;
x1F = (linspace(-1, 3.5, nPoints))';
x2F = linspace(-1, 3.5, nPoints);
FIM = zeros(nPoints);

maxFIM = 0;
maxFIMLocation = [];

for m1 = 1:nPoints
	for m2 = 1:nPoints
		expRange1	= norm([x1F(m1); x2F(m2)] - xHat);
		expRange2	= norm(sensor2 - xHat);
		if (expRange1 <= sensorMaxRange) && (expRange1 >= sensorMinRange)
			FIM(m2, m1) = (1 / rMeas) *  det( (1 / expRange1^2) * [...
				(xHat(1) - x1F(m1))^2						(xHat(1) - x1F(m1))*(xHat(2) - x2F(m2))			
				(xHat(1) - x1F(m1))*(xHat(2) - x2F(m2))		(xHat(2) - x2F(m2))^2 ] + ...
				(1 / expRange2^2) * [... 
				(xHat(1) - sensor2(1))^2					(xHat(1) - sensor2(1))*(xHat(2) - sensor2(2))			
				(xHat(1) - sensor2(1))*(xHat(2) - sensor2(2))		(xHat(2) - sensor2(2))^2
				]);
		end

		if FIM(m2, m1) >= maxFIM
			maxFIM = FIM(m2, m1);
			maxFIMLocation = [x1F(m1); x2F(m2)];
		end
	end
end
surf(ax, x1F, x2F, FIM, 'LineStyle','none');
colorbar; view(2);

zMax	= max(FIM(:));

xPlot	= linspace( (xHat(1) - (1-1E-6)*(sensorMinRange)), (xHat(1) + (1-1E-6)*(sensorMinRange)), 100 );
yPlot1	= xHat(2) + (sensorMinRange^2 - (xPlot - xHat(1)).^2).^0.5;
yPlot2	= xHat(2) - (sensorMinRange^2 - (xPlot - xHat(1)).^2).^0.5;
plot3(xPlot, yPlot1, zMax*ones(length(xPlot), 1), 'LineWidth', 2);
ax = gca; ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot3(xPlot, yPlot2, zMax*ones(length(xPlot), 1), 'LineWidth', 2);


xPlot	= linspace( (xHat(1) - (1-1E-6)*(sensorMaxRange)), (xHat(1) + (1-1E-6)*(sensorMaxRange)), 500 );
yPlot1	= xHat(2) + (sensorMaxRange^2 - (xPlot - xHat(1)).^2).^0.5;
yPlot2	= xHat(2) - (sensorMaxRange^2 - (xPlot - xHat(1)).^2).^0.5;
plot3(xPlot, yPlot1, zMax*ones(length(xPlot), 1), 'LineWidth', 2);
ax = gca; ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot3(xPlot, yPlot2, zMax*ones(length(xPlot), 1), 'LineWidth', 2);

plot3(maxFIMLocation(1), maxFIMLocation(2), zMax, '.', 'MarkerSize', 20)
plot3(sensor2(1), sensor2(2), zMax, 'o', 'MarkerSize', 10)

%(xHat(1) - x1F(m1))*(xHat(2) - x2F(m2))

return


%% Mutual information

maxMI = 0;
maxMILocation = [];

MI = zeros(nPoints);
for m1 = 1:nPoints
	for m2 = 1:nPoints
		expRange1	= norm([x1F(m1); x2F(m2)] - xHat);
		expRange2	= norm(sensor2 - xHat);

		C	= [...
			(1 / expRange1)*[(xHat(1) - x1F(m1))		(xHat(2) - x2F(m2)) ]; ...
			(1 / expRange2)*[(xHat(1) - sensor2(1))		(xHat(2) - sensor2(2)) ] ...
			];
		PZZ	= C * PXX * C' + rMeas*eye(2);
		PXZ	= PXX*C';
		PJoint = [PXX PXZ; PXZ' PZZ];

		if (expRange1 <= sensorMaxRange) && (expRange1 >= sensorMinRange)
			MI(m2, m1) = 0.5*log( det(PXX) * det(PZZ) / det(PJoint) );
		end

		if MI(m2, m1) >= maxMI
			maxMI = MI(m2, m1);
			maxMILocation = [x1F(m1); x2F(m2)];
		end
	end
end





zMax	= max(MI(:));


figure;
plot(xHat(1), xHat(2), 'x', 'MarkerSize', 10); hold on; axis equal
xPlot	= linspace( (xHat(1) - (1-1E-6)*sqrt(PXX(1,1))), (xHat(1) + (1-1E-6)*sqrt(PXX(1,1))), 100 );
yPlot1	= sqrt(PXX(2,2)) * (1 - (xPlot - xHat(1)).^2 / PXX(1,1) ).^0.5 + xHat(2);
yPlot2	= -sqrt(PXX(2,2)) * (1 - (xPlot - xHat(1)).^2 / PXX(1,1) ).^0.5 + xHat(2);

plot(xPlot, yPlot1, 'LineWidth', 2);
ax = gca; ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot(xPlot, yPlot2, 'LineWidth', 2);

make_nice_figures(gcf, gca, 14, [], '$x_1$', '$x_2$', 'Mutual Information', [], [], [-1 3.5], [-1 3.5])

surf(ax, x1F, x2F, MI, 'LineStyle','none');
colorbar; view(2);
plot3(maxMILocation(1), maxMILocation(2), zMax, '.', 'MarkerSize', 20)
plot3(sensor2(1), sensor2(2), zMax, 'o', 'MarkerSize', 10)


%% KL Divergence (2 sensors, 1 fixed)

KLDPosterior = zeros(nPoints);
maxKLD = 0;
maxKLDLocation = [];

for m1 = 1:nPoints
	for m2 = 1:nPoints
		expRange1	= norm([x1F(m1); x2F(m2)] - xHat);
		expRange2	= norm(sensor2 - xHat);

		C	= [...
			(1 / expRange1)*[(xHat(1) - x1F(m1))		(xHat(2) - x2F(m2)) ]; ...
			(1 / expRange2)*[(xHat(1) - sensor2(1))		(xHat(2) - sensor2(2)) ] ...
			];
		PZZ	= C * PXX * C' + rMeas*eye(2);
		PXZ	= PXX*C';
		L	= PXZ / PZZ;
		postP	= (eye(2) - L*C)*PXX;

		if (expRange1 <= sensorMaxRange) && (expRange1 >= sensorMinRange)
			KLDPosterior(m2, m1) = 0.5*( log( det(postP) / det(PXX) ) - 2 + ...
				trace((L'/postP)*L*PZZ) + trace(postP\PXX) );
		end

		if KLDPosterior(m2, m1) >= maxKLD
			maxKLD = KLDPosterior(m2, m1);
			maxKLDLocation = [x1F(m1); x2F(m2)];
		end
	end
end

maxMILocation
maxKLDLocation

zMax	= max(KLDPosterior(:));


figure;
plot(xHat(1), xHat(2), 'x', 'MarkerSize', 10); hold on; axis equal
xPlot	= linspace( (xHat(1) - (1-1E-6)*sqrt(PXX(1,1))), (xHat(1) + (1-1E-6)*sqrt(PXX(1,1))), 100 );
yPlot1	= sqrt(PXX(2,2)) * (1 - (xPlot - xHat(1)).^2 / PXX(1,1) ).^0.5 + xHat(2);
yPlot2	= -sqrt(PXX(2,2)) * (1 - (xPlot - xHat(1)).^2 / PXX(1,1) ).^0.5 + xHat(2);

plot(xPlot, yPlot1, 'LineWidth', 2);
ax = gca; ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot(xPlot, yPlot2, 'LineWidth', 2);

make_nice_figures(gcf, gca, 14, [], '$x_1$', '$x_2$', 'KL Divergence', [], [], [-1 3.5], [-1 3.5])

surf(ax, x1F, x2F, KLDPosterior, 'LineStyle','none');
colorbar; view(2);
plot3(maxKLDLocation(1), maxKLDLocation(2), zMax, '.', 'MarkerSize', 20)
plot3(sensor2(1), sensor2(2), zMax, 'o', 'MarkerSize', 10)

