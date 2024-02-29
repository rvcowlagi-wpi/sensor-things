clear variables; close all; clc

positionX0	= [10; 10];
timeStamps	= 0:0.01:25;

velocityV	= [5; 0];
positionX(1, :)	= positionX0(1) + timeStamps*velocityV(1);
positionX(2, :)	= positionX0(2) + timeStamps*velocityV(2);

rangeTrue	= (positionX(1, :).^2 + positionX(2, :).^2).^(0.5);

measNoiseRange = (0.5)*randn(1, length(timeStamps)); % Range error std. dev. is 0.5
rangeMeas	= rangeTrue + measNoiseRange;


save data_2D_RT_tracker.mat rangeMeas timeStamps