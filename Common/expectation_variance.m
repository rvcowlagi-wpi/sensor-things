clear variables; close all; clc

nSamples	= 1E6;
x	= 1 + rand(nSamples, 1);

mean(x)


y = 1 ./x;

mean(y)


