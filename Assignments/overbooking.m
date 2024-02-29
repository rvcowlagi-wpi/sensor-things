clear variables; close all; clc

n	= 100;
p	= 0.1;

%% Gamma Airlines
probOverbooking	= 0;
for k = 91:100
	probOverbooking = probOverbooking + ...
		nchoosek(n, k) * p^(n - k) * (1 - p)^k;
end
probOverbooking

%% Divided Airlines
n	= 200;
p	= 0.1;
lambda_ = n*p;
probOverbooking	= 0;
for k = 0:20
	probOverbooking = probOverbooking + ...
		(lambda_^k) * exp(-lambda_) / factorial(k);
end
probOverbooking
