function sim_multisensor_tracking()

%% Constants and parameters
close all; clc

simTimeStamps	= 0:1E-3:0.75;
nTimeStamps		= length(simTimeStamps);
dt_				= simTimeStamps(2) - simTimeStamps(1);

nSensors		= 2;
sensorLocations	= [rand(1, nSensors); zeros(1 ,nSensors)];
ntReconfig		= 10;
% sensorLocations is an array of two rows and Ns columns, where Ns is the
% number of sensors. Therefore, the k^th column in this array provides the
% coordinates of the k^th sensor. The initial sensor locations are spread
% out randomly along the y_1 axis. You are free to choose a different
% initialization.

%% Other problem data

stdDevMeasNoise = 0.03;
rangeMin		= 0.1;
rangeMax		= 1;
priorX0			= [1; 0];
priorP0			= 0.1*eye(2);

%% Main simulation loop
for m1 = 1:nTimeStamps

	t = simTimeStamps(m1);

	%------ Sensor reconfiguration after every ntReconfig time steps
	if ~mod(m1, ntReconfig)
		sensorLocations = sensor_configuration();							% <===== YOUR SENSOR RECONFIGURATION HERE; PLACEHOLDER FUNCTION BELOW
	end

	z	= asg5_sensor_emulator(t, sensorLocations);

	xHat= estimator();														% <===== YOUR ESTIMATION TECHNIQUE HERE; PLACEHOLDER FUNCTION BELOW
	
	%----- Store results to plot later
end



%%
	
	function newLocations_ = sensor_configuration()
		%**** YOUR SENSOR CONFIGURATION HERE *****
		%**** USE ANY INPUT ARGUMENTS TO THIS FUNCTION AS NECESSARY ****

		newLocations_ = sensorLocations;									% Placeholder
	end

	function xHat_ = estimator()
		%**** YOUR SENSOR CONFIGURATION HERE *****
		%**** USE ANY INPUT ARGUMENTS TO THIS FUNCTION AS NECESSARY ****
		%**** YOU MAY CHANGE THE OUTPUT ARGUMENTS TO INCLUDE ESTIMATION
		%		ERROR COVARIANCE, DEPENDING ON YOUR ESTIMATOR *****

		z;
		%** NOT IDEAL CODING PRACTICE BUT CONVENIENT, z IS AVAILABLE HERE
		% AS A "GLOBAL" VARIABLE WITHOUT HAVING TO PASS AS AN INPUT ARGUMENT
		xHat_ = [0; 0];
	end
end