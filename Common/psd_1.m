clear variables; close all; clc

fSample = 1e3;
t_		= 0:1/fSample:1-1/fSample;
% x = 0*cos(2*pi*100*t) + 5*randn(size(t));
x_		= cos(2*pi*10*t_) + sin(2*pi*100*t_);

% fs1 = 10;
% t1 = 0:1/fs1:1;

% x1 = randn(size(t1));
% x_ = interp1(t1,x1,t_,"nearest");

N		= length(x_);
xDFT	= fft(x_);
xDFT	= xDFT(1:N/2+1);
psdx	= (1/(fSample*N)) * abs(xDFT).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fSample/length(x_):fSample/2;

plot(freq,pow2db(psdx))
grid on
title("Periodogram Using FFT")
xlabel("Frequency (Hz)")
ylabel("Power/Frequency (dB/Hz)")
