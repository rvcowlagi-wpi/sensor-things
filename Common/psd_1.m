function psd_1()

close all; clc

fSample = 1e3;
time_		= 0:1/fSample:1-1/fSample;
% xSignal_	= cos(2*pi*10*t_) + sin(2*pi*100*t_)
xSignal_	= randn(size(time_));
[~, y1]		= ode45(@(t,x) simpleLTI(t,x, xSignal_), time_, 0);	
ySignal_	= y1(:, 1)'; 

% fs1 = 10;
% t1 = 0:1/fs1:1;

% x1 = randn(size(t1));
% x_ = interp1(t1,x1,t_,"nearest");

N		= length(ySignal_);
xDFT	= fft(ySignal_);
xDFT	= xDFT(1:N/2+1);
psdx	= (1/(fSample*N)) * abs(xDFT).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fSample/length(ySignal_):fSample/2;

% fig0	= figure;
% plot(time_, xSignal_, 'LineWidth', 3); grid on;
% ax0		= gca;
% make_nice_figures(fig0, ax0, 18, [], ...
% 	'Time (s)', 'Signal $x(t)$', 'Signal', ...
% 	[], 'psd_noise_signal.png', [], [-3 3])
% 
% fig1	= figure;
% plot(freq, pow2db(psdx), 'LineWidth', 3); grid on;
% ax1		= gca;
% make_nice_figures(fig1, ax1, 18, [], ...
% 	'Frequency (Hz)', 'Power/Frequency (dB/Hz)', 'Periodogram', ...
% 	[], 'psd_noise_pdg.png', [], [-350 0])

fig0	= figure;
plot(time_, ySignal_, 'LineWidth', 3); grid on;
ax0		= gca;
make_nice_figures(fig0, ax0, 18, [], ...
	'Time (s)', 'Signal $y(t)$', 'Signal', ...
	[], 'psd_preshaped2_noise_signal.png', [], [])

fig1	= figure;
plot(freq, pow2db(psdx), 'LineWidth', 3); grid on;
ax1		= gca;
make_nice_figures(fig1, ax1, 18, [], ...
	'Frequency (Hz)', 'Power/Frequency (dB/Hz)', 'Periodogram', ...
	[], 'psd_preshaped2_noise_pdg.png', [], [])


	function xDot_ = simpleLTI(t_, x_, u_)
		tIndex	= max(1, min(round(t_*fSample), length(time_)));
		ut_		= u_(tIndex);

% 		xDot_	= [0 1; -4 -5]*x_ + [0; 1]*ut_;
		xDot_	= -2*x_ + ut_; % Brown noise
	end

end