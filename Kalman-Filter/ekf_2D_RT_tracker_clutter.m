%% Extended Kalman Filter implementation for 2D Range-Bearing Tracking
% Copy-pasted from 1D navigation; fix comments later

function ekf_2D_RT_tracker_clutter()

close all; clc

%%
load data_2D_RT_tracker_clutter.mat zRange zRange2 zRange3 ...
	zBear zBear2 zBear3 timeStamps rangeTrue bearTrue

nStates = 4;	% States are r, rDot, theta, thetaDot
nMeas	= 2;	% Measurements are r, theta

V		= 50;

%% Error Covariances
Q	= [(1E-2)^2 0; 0 (1E-2*pi/180)^2];
R	= [0.3^2 0; 0 (2*pi/180)^2];

%% Time Step and Interval of Interest
%%
% We have uniformly spaced time stamps, so calculating the difference
% between any two successive instants gives us $dt.$
dt	= timeStamps(2) - timeStamps(1);
%%
% Note the number of time stamps, which will be the number of iterations
% for which the Kalman filter will run.
nTimeStamps = length(timeStamps);

%% Predictive Model and Measurement Model
G	= [ [0; dt] zeros(2, 1); zeros(2, 1) [0; dt] ];
C	= [1 0 0 0; 0 0 1 0];

%% Initialization
xHat	= [zRange(1); V; zBear(1); -V/zRange(1)];	
P		= diag([R(1,1), 0.1, R(2,2), 0.1]);
% 
% [PSteadyState,~,~]	= idare(A', C', (G*Q*G'), R, [], []);
% PSteadyStateStore	= reshape(PSteadyState, nStates^2, 1);
% LSteadyState		= PSteadyState * C' / (C*PSteadyState*C' + R);
%% 
storeXHat	= zeros(nStates, nTimeStamps);
storeP		= zeros(nStates^2, nTimeStamps);
storePTrace	= zeros(1, nTimeStamps);
storeInnov	= zeros(nMeas, nTimeStamps);
storeInnovCovar	= zeros(nMeas^2, nTimeStamps);
storeNormInnov	= zeros(1, nTimeStamps);

%%
storeXHat(:, 1)		= xHat;
storePTrace(:, 1)	= trace(P);
storeP(:, 1)		= reshape(P, nStates^2, 1); 

%% Run Kalman Filter Iterations
for m1 = 1:(nTimeStamps-1)

	%----- Prediction equations to get a preliminary estimate and error covariance.
	u			= 0;
	A			= eye(nStates) + jacobianA(xHat)*dt;

% 	xHatMinus	= xHat + polar_kinematics_2D([], xHat, V)*dt;				% First-order Euler approximation
% 	[~, xSim]	= ode45(@(t,x) polar_kinematics_2D(t,x, V), [0 dt], xHat);	% RK4
% 	xHatMinus	= xSim(end, :)';
	xHatMinus	= one_step_update(xHat, [0; 0]);
	PMinus		= A*P*A' + G*Q*G';

	%----- Compute Kalman gain. Note the use of |/| command instead of inverse.
	L			= (PMinus * C') / (C * PMinus * C' + R );
	
	%----- Measurement update
	z			= [zRange(m1 + 1); zBear(m1 + 1)];
	thisInnov	= z - C*xHatMinus;	

	thisInnovCovar		= C * PMinus * C' + R;
	thisNormalizedInnov = ( thisInnov' / thisInnovCovar ) * thisInnov;

	if m1 == 25 
		figure;
		cla;
		plot(rangeTrue(m1)*cos(bearTrue(m1)), rangeTrue(m1)*sin(bearTrue(m1)), ...
			'.', 'MarkerSize', 20); hold on;
		
		plot(zRange(m1 + 1)*cos(zBear(m1 + 1)), zRange(m1 + 1)*sin(zBear(m1 + 1)), ...
			'x', 'MarkerSize', 10); 
		ax = gca; ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
		hold on;

		plot(zRange2(m1 + 1)*cos(zBear2(m1 + 1)), zRange2(m1 + 1)*sin(zBear2(m1 + 1)), ...
			'x', 'MarkerSize', 10); ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
		plot(zRange3(m1 + 1)*cos(zBear3(m1 + 1)), zRange3(m1 + 1)*sin(zBear3(m1 + 1)), ...
			'x', 'MarkerSize', 10)
	
		plot(storeXHat(1, 1:m1).*cos(storeXHat(3, 1:m1)), ...
			storeXHat(1, 1:m1).*sin(storeXHat(3, 1:m1)), ...
			'.', 'MarkerSize', 30)
		make_nice_figures(gcf, gca, 18, [], '$y_1$', '$y_2$', ...
			'Range Rate', [],[],[0 80],[0 20]);
		drawnow()
% 		exportgraphics(gca,'mht-example1.png','Resolution',150)
		pause(0.1)

		make_nice_figures(gcf, gca, 18, [], '$y_1$', '$y_2$', ...
			'Range Rate', [],[],[21 25], [9 11]);
% 		exportgraphics(gca,'mht-example2.png','Resolution',150)


		z1	= [zRange(m1 + 1); zBear(m1 + 1)]
		z2	= [zRange2(m1 + 1); zBear2(m1 + 1)]
		z3	= [zRange3(m1 + 1); zBear3(m1 + 1)]


		thisInnov1	= z1 - C*xHatMinus
		thisInnov2	= z2 - C*xHatMinus
		thisInnov3	= z3 - C*xHatMinus

		( thisInnov1' / thisInnovCovar ) * thisInnov1
		( thisInnov2' / thisInnovCovar ) * thisInnov2
		( thisInnov3' / thisInnovCovar ) * thisInnov3

		xHat1		= xHatMinus + L*(thisInnov1)
		xHat2		= xHatMinus + L*(thisInnov2)

		plot( xHat1(1)*cos(xHat1(3)), xHat1(1)*sin(xHat1(3)), '.', 'MarkerSize', 30)
		text( xHat1(1)*cos(xHat1(3)), xHat1(1)*sin(xHat1(3)), '$\hat{x}_1$', 'Interpreter','latex')
		plot( xHat2(1)*cos(xHat2(3)), xHat2(1)*sin(xHat2(3)), '.', 'MarkerSize', 30)
		text( xHat2(1)*cos(xHat2(3)), xHat2(1)*sin(xHat2(3)), '$\hat{x}_2$', 'Interpreter','latex')
% 		exportgraphics(gca,'mht-example3.png','Resolution',150)

		lambda1		= ( thisInnov1' / thisInnovCovar ) * thisInnov1
		lambda2		= ( thisInnov2' / thisInnovCovar ) * thisInnov2
% 		xHat3		= xHatMinus + L*(thisInnov3);

	end

	if m1 == 26

		xHatMinus1	= one_step_update(xHat1, [0; 0]);
		xHatMinus2	= one_step_update(xHat2, [0; 0]);

		z1	= [zRange(m1 + 1); zBear(m1 + 1)]
		z2	= [zRange2(m1 + 1); zBear2(m1 + 1)]
		z3	= [zRange3(m1 + 1); zBear3(m1 + 1)]


		thisInnov11	= z1 - C*xHatMinus1
		thisInnov21	= z2 - C*xHatMinus1
		thisInnov31	= z3 - C*xHatMinus1
		thisInnov12	= z1 - C*xHatMinus2
		thisInnov22	= z2 - C*xHatMinus2
		thisInnov32	= z3 - C*xHatMinus2

		( thisInnov11' / thisInnovCovar ) * thisInnov11
		( thisInnov21' / thisInnovCovar ) * thisInnov21
		( thisInnov31' / thisInnovCovar ) * thisInnov31
		( thisInnov12' / thisInnovCovar ) * thisInnov12
		( thisInnov22' / thisInnovCovar ) * thisInnov22
		( thisInnov32' / thisInnovCovar ) * thisInnov32

		xHat1		= xHatMinus1 + L*(thisInnov11)
		xHat2		= xHatMinus1 + L*(thisInnov21)
		xHat3		= xHatMinus2 + L*(thisInnov12)
		xHat4		= xHatMinus2 + L*(thisInnov22)

		plot( xHat1(1)*cos(xHat1(3)), xHat1(1)*sin(xHat1(3)), '.', 'MarkerSize', 30)
		text( xHat1(1)*cos(xHat1(3)), xHat1(1)*sin(xHat1(3)), '$\hat{x}_1$', 'Interpreter','latex')

		plot( xHat2(1)*cos(xHat2(3)), xHat2(1)*sin(xHat2(3)), '.', 'MarkerSize', 30)
		text( xHat2(1)*cos(xHat2(3)), xHat2(1)*sin(xHat2(3)), '$\hat{x}_2$', 'Interpreter','latex')

		plot( xHat3(1)*cos(xHat3(3)), xHat3(1)*sin(xHat3(3)), '.', 'MarkerSize', 30)
		text( xHat3(1)*cos(xHat3(3)), xHat3(1)*sin(xHat3(3)), '$\hat{x}_3$', 'Interpreter','latex')


		plot( xHat4(1)*cos(xHat4(3)), xHat4(1)*sin(xHat4(3)), '.', 'MarkerSize', 30)
		text( xHat4(1)*cos(xHat4(3)), xHat4(1)*sin(xHat4(3)), '$\hat{x}_4$', 'Interpreter','latex')

% 		exportgraphics(gca,'mht-example4.png','Resolution',150)

		lambda11	= lambda1 + ( thisInnov11' / thisInnovCovar ) * thisInnov11
		lambda12	= lambda1 + ( thisInnov21' / thisInnovCovar ) * thisInnov21
		lambda21	= lambda2 + ( thisInnov12' / thisInnovCovar ) * thisInnov12
		lambda22	= lambda2 + ( thisInnov22' / thisInnovCovar ) * thisInnov22


	end

	xHat		= xHatMinus + L*(thisInnov);
	P			= (eye(nStates) - L*C)*PMinus;
	
	%----- Store results
	storeXHat(:, m1+1)		= xHat;
	storeP(:, m1+1)			= reshape(P, nStates^2, 1);
	storePTrace(:, m1+1)	= trace(P);
	storeInnov(:, m1+1)		= thisInnov;
	storeInnovCovar(:, m1+1)= reshape(thisInnovCovar, nMeas^2, 1);
	storeNormInnov(:, m1+1)	= thisNormalizedInnov;
end

%% Innovation autocorrelation
nTimeStamps		= length(timeStamps);
nMeasDim		= 2;

innovAutoCorr	= zeros(1, nTimeStamps);
shiftedInnov	= zeros(nMeasDim, nTimeStamps);
twoSigmaBd		= zeros(1, nTimeStamps);
for m1 = 1:nTimeStamps - 10
	shiftedInnov(:,m1+1:end)= storeInnov(:, 2:(nTimeStamps - m1 + 1));
	innovAutoCorr(m1)		= (1 / (nTimeStamps - m1) ) * (...
			sum( shiftedInnov(1, m1+1:end).*storeInnov(1, m1+1:end) ) + ...
			sum( shiftedInnov(2, m1+1:end).*storeInnov(2, m1+1:end) ) );
	twoSigmaBd(m1)			= 2/sqrt(nTimeStamps - m1); 
end

%% Figures of merit of the filter
% format shortE

mseeRange	= mean( (storeXHat(1,:) - rangeTrue).^2 );
mseeBrng	= mean( (storeXHat(3,:) - bearTrue).^2 );

fprintf('----- Mean of innovations autocorrelation (closer to zero is better): \n'); disp( mean(innovAutoCorr(2:end)) )
fprintf('----- Mean square estimation error in range (closer to zero is better): \n'); disp( mseeRange )
fprintf('----- Mean square estimation error in bearing (closer to zero is better): \n'); disp( mseeBrng )


%% Plot Results
figure;
for m1 = 1:nTimeStamps
	cla;
	plot(rangeTrue(m1)*cos(bearTrue(m1)), rangeTrue(m1)*sin(bearTrue(m1)), ...
		'.', 'MarkerSize', 20); hold on;
	ax = gca;
	plot(zRange(m1)*cos(zBear(m1)), zRange(m1)*sin(zBear(m1)), ...
		'x', 'MarkerSize', 10); ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
	plot(zRange2(m1)*cos(zBear2(m1)), zRange2(m1)*sin(zBear2(m1)), ...
		'x', 'MarkerSize', 10); ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
	plot(zRange3(m1)*cos(zBear3(m1)), zRange3(m1)*sin(zBear3(m1)), ...
		'x', 'MarkerSize', 10)

	plot(storeXHat(1, 1:m1).*cos(storeXHat(3, 1:m1)), ...
		storeXHat(1, 1:m1).*sin(storeXHat(3, 1:m1)), ...
		'.', 'MarkerSize', 10)
	make_nice_figures(gcf, gca, 18, [], '$y_1$', '$y_2$', ...
		'Range Rate', [],[],[0 80],[0 20]);
	drawnow()
	pause(0.1)
	
end

% fig1 = figure; 
% subplot(221); hold on; plot(timeStamps, storeXHat(1,:), 'LineWidth', 2);
% hold on;
% plot(timeStamps, rangeTrue, 'LineWidth',  2)
% make_nice_figures(gcf, gca, 18, [], 'Time (h)', ...
% 	'Range $\hat{x}_1 = r$ (km)', 'Range and Rate', [],[],[],[]);
% 
% subplot(222); hold on; plot(timeStamps, storeXHat(2,:), 'LineWidth', 2)
% make_nice_figures(gcf, gca, 18, [], 'Time (h)', ...
% 	'Range rate $\hat{x}_2 = \dot{r}$ (km/hr)', 'Range and Rate', [],[],[],[]);
% 
% subplot(223); hold on; plot(timeStamps, storeXHat(3,:)*180/pi, 'LineWidth', 2);
% hold on;
% plot(timeStamps, bearTrue*180/pi, 'LineWidth',  2)
% make_nice_figures(gcf, gca, 18, [], 'Time (h)', ...
% 	'Range $\hat{x}_3 = \theta$ (deg)', 'Bearing and Rate', [],[],[],[]);
% 
% subplot(224); hold on; plot(timeStamps, storeXHat(4,:)*180/pi, 'LineWidth', 2)
% make_nice_figures(gcf, gca, 18, [], 'Time (h)', ...
% 	'Range rate $\hat{x}_4 = \dot{\theta}$ (deg/hr)', 'Bearing and Rate', [],[],[],[]);
% 
% 
% 
% 
% fig2 = figure;
% plot(timeStamps, storePTrace, 'LineWidth', 2);
% make_nice_figures(gcf, gca, 18, [], 'Time (h)', 'tr$(P)$', 'Trace', [],[],[],[]);
% 
% fig3 = figure;
% plot(timeStamps, storeP, 'LineWidth', 2); hold on;
% make_nice_figures(gcf, gca, 18, [], 'Time (h)', '$p_{ij}$', 'Trace', [],[],[],[]);
% 
% ylabel('$p_{ij}$', 'Interpreter', 'latex', 'FontSize', 12);
% 
% 
% figure;
% subplot(121)
% plot(timeStamps, storeInnov(1, :), 'LineWidth', 2); hold on;
% ax = gca; ax.ColorOrderIndex = 1;
% plot(timeStamps, 2*(storeInnovCovar(1,:).^0.5), ...
% 	'LineWidth', 2, 'LineStyle', '--'); ax.ColorOrderIndex = 1;
% plot(timeStamps, -2*(storeInnovCovar(1,:).^0.5), ...
% 	'LineWidth', 2, 'LineStyle', '--')
% 
% subplot(122)
% plot(timeStamps, storeInnov(2, :), 'LineWidth', 2); hold on;
% ax = gca; ax.ColorOrderIndex = 2;
% plot(timeStamps, 2*(storeInnovCovar(4, :).^0.5), ...
% 	'LineWidth', 2, 'LineStyle', '--'); ax.ColorOrderIndex = 2;
% plot(timeStamps, 2*(storeInnovCovar(4, :).^0.5), ...
% 	'LineWidth', 2, 'LineStyle', '--')
% 
% title('Innovations sequence', 'Interpreter', 'latex', 'FontSize', 14)
% xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 12);
% ylabel('Innovation', 'Interpreter', 'latex', 'FontSize', 12);
% 
% 
% figure;
% plot(timeStamps, innovAutoCorr, 'LineWidth', 2); hold on;
% plot(timeStamps, (2/sqrt(nTimeStamps))*ones(size(timeStamps)), 'LineWidth', 2, 'LineStyle', '--' );
% ax = gca; ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
% plot(timeStamps, (-2/sqrt(nTimeStamps))*ones(size(timeStamps)), 'LineWidth', 2, 'LineStyle', '--' )
% title('Innovations autocorrelation', 'Interpreter', 'latex', 'FontSize', 14)
% xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 12);
% ylabel('Innovation', 'Interpreter', 'latex', 'FontSize', 12);
% 
% 
% figure;
% histogram(storeNormInnov, 'Normalization','pdf', 'NumBins', 10); hold on;
% xPlot = linspace(0, 15, 100);
% plot( xPlot, chi2pdf(xPlot, nMeas), 'LineWidth', 3);
% make_nice_figures(gcf, gca, 14, [], 'Normalized innovation', 'PDF', 'Normalized Innovation', [],[],[],[]);
% 
% figure;
% plot(xPlot, chi2cdf(xPlot, nMeas), 'LineWidth', 3); grid on;
% make_nice_figures(gcf, gca, 14, [], 'Normalized innovation', 'CDF', 'Chi-square CDF', [],[],[],[]);


%%
	function xDot = system_dynamics_continuous(t_, x_, w_)
		xDot	= zeros(4, 1);
	
	% 	xDot(1) = x_(2);
	% 	xDot(3)	= x_(4);
		xDot(1)	= V*cos(x_(3));
		xDot(3)	= -V*sin(x_(3)) / x_(1);
		xDot(2)	= -V*x_(4)*sin(x_(3)) + w_(1);
		xDot(4)	= V*x_(2)*sin(x_(3)) / (x_(1)^2) - V*x_(4)*cos(x_(3)) / x_(1) + w_(2);
	end
	
	function jacA = jacobianA(x_)
		jacA = zeros(4);
		jacA(1, :)	= [0 1 0 0];
		jacA(3, :)	= [0 0 0 1];
		jacA(2, :)	= V*[0 0 -x_(4)*cos(x_(3)) -sin(x_(3))];
		jacA(4, :)	= V*[ ...
			-2*x_(2)*sin(x_(3))/(x_(1)^3) + x_(4)*cos(x_(3))/(x_(1)^2) ...
			sin(x_(3))/(x_(1)^2) ...
			x_(2)*cos(x_(3))/(x_(1)^2) + x_(4)*sin(x_(3))/x_(1) ...
			-cos(x_(3)) / x_(1)];
	end

	function xk_ = one_step_update(xkMinus1_, wkMinus1_)
		xk_ = xkMinus1_ + ...
			system_dynamics_continuous([], xkMinus1_, wkMinus1_)*dt;		% First-order Euler approximation

		% Alternative to Euler approximation: RK4 (need to fix code to
		% accommodate w)
		% [~, xSim]	= ode45(@polar_kinematics_2D(t,x), [0 dt], xkMinus1_);
		% xk_	= xSim(end, :)';
	end

end