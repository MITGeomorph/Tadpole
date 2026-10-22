function Snapshot(p,g)

pfinal = SetUpPlot(p,g);

% hide the axes toolbars, which MATLAB warns about when printing a figure
ax = findall(pfinal.fighandle, 'Type', 'axes');
for k = 1:numel(ax)
    if isprop(ax(k), 'Toolbar') && ~isempty(ax(k).Toolbar), ax(k).Toolbar.Visible = 'off'; end
end

print(pfinal.fighandle, '-dpng', p.runname)

close(pfinal.fighandle)
