function range_ = asg5_sensor_emulator(t_, sensorLocations)


	load assignment5_problem2.mat timeStamps x

	x		= (x).^(1/2) - 25;
	tIndex	= find(t_ >= timeStamps, 1, 'last');

	rangeMin	= 0.1;
	rangeMax	= 1;

	xt		= x(:, tIndex);

	nSensors= size(sensorLocations, 2);
	range_	= zeros(nSensors, 1);

	for m1 = 1:nSensors
		thisRange	= norm(xt - sensorLocations(:, m1));
		if thisRange < rangeMin
			range_(m1)	= rangeMin + 0.03*randn;
		elseif thisRange > rangeMax
			range_(m1)	= rangeMax + 0.03*randn;
		else 
			range_(m1)	= thisRange + 0.03*randn;
		end
	end

end
