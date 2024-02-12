function BREx_latency()

close all; clc

a = 5;
z = 0.05; %a*rand

tSS		= linspace(0, a, 300);
priorT = (1/a)*ones(1, 300);
MAPEstimate = 1/z;


fontSize_ = 40;
figure('units', 'normalized', 'OuterPosition', [0.05 0.05 0.8 0.9]);

plot(tSS, priorT, tSS, posterior_TGivenZ(tSS, z), 'LineWidth', 3); hold on;
plot(z*ones(1, 100), linspace(0, max(posterior_TGivenZ(tSS, z)), 100), ...
	'LineStyle','--', 'LineWidth', 3 )
plot(MAPEstimate*ones(1, 100), linspace(0, ...
	max(posterior_TGivenZ(tSS, z)), 100), 'LineStyle','--', 'LineWidth', 3 )

legend('Prior', 'Posterior', '$z$', 'MAP', 'interpreter', 'latex')

xlabel('$\theta$', 'FontName', 'Times New Roman', ...
	'FontSize', fontSize_, 'FontWeight', 'bold', 'interpreter', 'latex');
ylabel('$f_{\Theta | Z}$', 'FontName', 'Times New Roman', ...
	'FontSize', fontSize_, 'FontWeight', 'bold', 'interpreter', 'latex');
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = fontSize_;


%---- Sanity check
integral(@(t) posterior_TGivenZ(t, z), 0, a)

	function pdf_ = posterior_TGivenZ(t_, z_)
		pdf_ = t_ .* exp(-t_.*z_) .* (z_.^2) ./  (1 - exp(-a*z_).* (1 + a*z_) );
	end

end