function handles = EnableControls(handles,Option)
switch Option
    case {'IRLab2','IRLab1'}
        handles.BkgSubTick.Enable           = 'On';
        handles.mintimeBkg.Enable           = 'On';
        handles.toBkgtext.Enable            = 'On';
        handles.maxtimeBkg.Enable           = 'On';
        handles.unitsBkgsub.Enable          = 'On'; 
        handles.PlotBkg.Enable              = 'On';
        handles.PlotNoise.Enable            = 'On';
        handles.make2Dplot.Enable           = 'On';
        handles.make3Dplot.Enable           = 'On';
        handles.plotKin.Enable              = 'On';
        handles.plotDeltaAbs.Enable         = 'On';
        handles.AddReplace1Dplot.Enable     = 'On';
        handles.InteractiveModeTick.Enable  = 'On';
        handles.doSVD.Enable                = 'On';
        handles.SVD_fit.Enable              = 'On';
        handles.NormKin.Enable              = 'On';
        handles.AnisotropyCalc.Enable      	= 'On';
        handles.SaveButton.Enable           = 'On';
        handles.Plot_Dependency.Enable      = 'On';
        handles.Fit.Enable                  = 'On';
        handles.IncludeSteadyState.Enable   = 'On';
        handles.tWLShift.Enable             = 'On';
        handles.nSVD.Enable                 = 'On';
        handles.nSVD.Visible                = 'On';
        handles.SaveTraces.Enable           = 'On';
        handles.PlotSumDiff.Enable          = 'On';
        handles.OffsetSubtraction.Enable    = 'On';
        handles.KineticsPerScan.Enable      = 'On';
        handles.SpectraPerScan.Enable       = 'On';
        handles.BinScans.Visible            = 'On';
        handles.BinScans.Enable             = 'On';
        handles.text35.Visible              = 'On';
        handles.text36.Visible              = 'On';
        handles.text35.Enable               = 'On';
        handles.text36.Enable               = 'On';
    case 'TRES'
        handles.BkgSubTick.Enable           = 'On';
        handles.mintimeBkg.Enable           = 'On';
        handles.toBkgtext.Enable            = 'On';
        handles.maxtimeBkg.Enable           = 'On';
        handles.unitsBkgsub.Enable          = 'On'; 
        handles.PlotBkg.Enable              = 'On';
        handles.PlotNoise.Enable            = 'Off';
        handles.make2Dplot.Enable           = 'On';
        handles.make3Dplot.Enable           = 'On';
        handles.plotKin.Enable              = 'On';
        handles.plotDeltaAbs.Enable         = 'On';
        handles.AddReplace1Dplot.Enable     = 'On';
        handles.InteractiveModeTick.Enable  = 'On';
        handles.doSVD.Enable                = 'On';
        handles.SVD_fit.Enable              = 'On';
        handles.KineticsPerScan.Enable      = 'Off';
        handles.SpectraPerScan.Enable       = 'Off';
        handles.NormKin.Enable              = 'On';
        handles.AnisotropyCalc.Enable  		= 'On';
        handles.SaveButton.Enable           = 'On';
        handles.Plot_dependency.Enable      = 'Off';
        handles.Fit.Enable                  = 'On';
        handles.IncludeSteadyState.Enable 	= 'Off';
        handles.tWLShift.Enable             = 'Off';
        handles.nSVD.Enable                 = 'On';
        handles.nSVD.Visible                = 'On';
        handles.SaveTraces.Enable           = 'On';
        handles.PlotSumDiff.Enable          = 'On';
        handles.OffsetSubtraction.Enable    = 'On';
        handles.BinScans.Visible            = 'Off';
        handles.BinScans.Enable             = 'Off';
        handles.text35.Visible              = 'Off';
        handles.text36.Visible              = 'Off';
        handles.text35.Enable               = 'Off';
        handles.text36.Enable               = 'Off';
end
