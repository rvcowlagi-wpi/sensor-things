function BREx_localization()

nSensors = 4;
sensorLocations = zeros(2,nSensors);

sensorLocations(:, 1) = [1; 1];
sensorLocations(:, 2) = [-1; 1];
sensorLocations(:, 3) = [-1; -1];
sensorLocations(:, 4) = [1; -1];

xTarget = -1 + rand(2, 1);

distances_ = zeros(1, nSensors);
for k = 1:nSensors
	distances_(k) = norm(sensorLocations(:, k) - xTarget);
end

detections_ = zeros(1, nSensors);
for k = 1:nSensors
	p_				= exp( -distances_(k) );
	detections_(k)	= binornd(1, p_);
end

	function pdf_ = joint_XZ(x_, z_)
		joint_XZ = prior_X(x_)
	end


end