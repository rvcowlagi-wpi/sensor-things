function make_nice_figures(fig_, ax_, fontSize_, axTitle_, ...
	xLabel_, yLabel_, figTitle_, figPosition_, exportAs_, xLim_, yLim_)

fig_.Units		= 'normalized';
if numel(figTitle_)
	fig_.Name		= figTitle_;
end
if numel(figPosition_)
	fig_.Position	= figPosition_;
else
	fig_.Position	= [0.1 0.1 0.5*[1 1]];
end

ax_.FontSize	= fontSize_; 
ax_.FontName	= 'Times New Roman';
ax_.TickLabelInterpreter = 'latex';
if numel(axTitle_)
	ax_.Title.String= axTitle_;
end
if numel(xLabel_)
	ax_.XLabel.String = xLabel_;
	ax_.XLabel.Interpreter = 'latex';
end
if numel(yLabel_)
	ax_.YLabel.String = yLabel_;
	ax_.YLabel.Interpreter = 'latex';
end
if numel(xLim_)
	ax_.XLim = xLim_;
end
if numel(yLim_)
	ax_.YLim = yLim_;
end


if numel(exportAs_)
	exportgraphics(fig_, exportAs_, 'Resolution',150);
end