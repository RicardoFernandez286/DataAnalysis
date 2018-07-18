function EnableControls_2DIR(handles)

% Processing controls
handles.InteractiveModeTick.Enable  = 'On';
handles.BkgSubTick.Enable           = 'On';
handles.ApplyChanges.Enable         = 'On';
handles.Probe_Calibration.Enable    = 'On';
handles.PumpCorrection_tick.Enable  = 'On';
handles.make2Dplot.Enable           = 'On';
handles.make3Dplot.Enable           = 'On';
handles.plot_Kinetics.Enable        = 'On';
handles.plot_Slices.Enable          = 'On';
handles.Normalise.Enable			= 'On';

% Preview window controls
handles.ShowTimeDomain.Enable		= 'On';
handles.PixelNumber.Enable			= 'On';
handles.ShowPhasing.Enable			= 'On';
handles.ShowProbeCalibration.Enable	= 'On';
handles.SaveProbeCal.Enable         = 'On';

% Phasing, zeropadding and apodisation controls
handles.phase_method.Enable			= 'On';
handles.phase_Npoints.Enable	    = 'On';
handles.apodise_method.Enable	    = 'On';
handles.zeropad_tick.Enable		    = 'On';
handles.zeropad_factor.Enable	    = 'On';
handles.zeropad_next2k.Enable	    = 'On';
handles.phase_coeffs.Enable			= 'Inactive';
handles.text31.Enable               = 'On';
handles.text32.Enable               = 'On';
handles.text34.Enable               = 'On';
handles.text43.Enable               = 'On';

% Analysis controls
handles.IntegralDynamics.Enable		= 'On';
handles.ShiftT2.Enable              = 'On';
handles.SubtractSpectra.Enable     = 'On';
handles.SpectralDiffusion.Enable    = 'On';
handles.Fit.Enable                  = 'On';
handles.PowerDependence.Enable      = 'On';
