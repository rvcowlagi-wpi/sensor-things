clear variables; close all; clc

x1True	= 0.1;
x2True	= 1.5;

r		= 0.01;

t_		= (0:1:4)';

z_		= [ones(5, 1) t_]*[x1True; x2True] + sqrt(r)*randn(5, 1);

%---- With all "batch" measurements
xBar	= [0; 0];
PXX		= [0.1 0; 0 0.1];

C	= [ones(5, 1) t_];
PZ	= C*PXX*C' + r*eye(5);
PXZ	= (C*PXX)';

xHatBatch	= xBar + PXZ*( PZ \ (z_ - C*xBar) )
PX			= PXX - (PXZ / PZ)*(PXZ')



%---- With two "batch" measurements
xBar	= [0; 0];
PXX		= [0.1 0; 0 0.1];

C = [ [1; 1] t_(1:2)];
PZ	= C*PXX*C' + r*eye(2);
PXZ	= (C*PXX)';

xHat2Batch	= xBar + PXZ*( PZ \ (z_(1:2) - C*xBar) )
PX			= PXX - (PXZ / PZ)*(PXZ')


%----- Recursive
xBar	= [0; 0];
PXX		= [0.1 0; 0 0.1];
xHat	= xBar;

for m = 1:length(t_)
	C = [1 t_(m)];
	PZ	= C*PXX*C' + r;
	PXZ	= (C*PXX)';

	xHat	= xHat + PXZ*( PZ \ (z_(m) - C*xHat) )
	PXX		= PXX - (PXZ / PZ)*(PXZ')

end
