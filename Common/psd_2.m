clear variables; close all; clc

dt	= 1;
t_	= 0:dt:10;

x_	= sin(5*t_);
periodogram(x_);

% [px, omga] = periodogram(x_);

plot(omga, 10*log10(px))