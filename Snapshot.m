function Snapshot(p,g)

pfinal = SetUpPlot(p,g);

print(pfinal.fighandle, '-dpng', p.runname)

close(pfinal.fighandle)
