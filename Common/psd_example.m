clear variables; close all; clc

fs = 1e4;
t = 0:1/fs:1-1/fs;
% x = 0*cos(2*pi*100*t) + 5*randn(size(t));
x = cos(2*pi*1000*t) + sin(2*pi*2000*t);

fs1 = 10;
t1 = 0:1/fs1:1;

x1 = randn(size(t1));
x = interp1(t1,x1,t,"nearest");

N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/length(x):fs/2;

plot(freq,pow2db(psdx))
grid on
title("Periodogram Using FFT")
xlabel("Frequency (Hz)")
ylabel("Power/Frequency (dB/Hz)")