function handles = UpdateEdits_2DIR(handles,Action)
minwl       = num2str(min(handles.ProbeAxis));
maxwl       = num2str(max(handles.ProbeAxis));

m           = handles.Population_delay.Value;
k           = 1; % Show the first datastate

Ncontours   = '40';
plot_cRange = '100';
plot_percentwhites = '5';

minindex    = findClosestId2Val(handles.PumpAxis{m,k},min(handles.ProbeAxis));
maxindex    = findClosestId2Val(handles.PumpAxis{m,k},max(handles.ProbeAxis));

maxabs      = max(max(handles.PROC_2D_DATA{m,k}(minindex:maxindex,:)));
minabs      = min(min(handles.PROC_2D_DATA{m,k}(minindex:maxindex,:)));
Zminmax     = [num2str(minabs,'%.4g') ', ' num2str(maxabs,'%.4g')];

% Show resolution statistics
resW1       = handles.PumpAxis{m,k}(2)-handles.PumpAxis{m,k}(1);
resW3       = handles.ProbeAxis(2)-handles.ProbeAxis(1);
handles.omega1_resolution_text.String = num2str(resW1);
handles.omega3_resolution_text.String = num2str(resW3);

% Display bin zero, phase and FT size
handles.binzero_text.String     = num2str(handles.binzero{m,k});
handles.phase_text.String       = num2str(handles.ZP_phase{m,k}(handles.binspecmax(m,k)));
handles.FTsize_text.String      = num2str(length(handles.PumpAxis{m,k}));

% Show number of scans
handles.Nscans_number.String = num2str(handles.Nscans);

% Show the phase coefficients
handles.phase_coeffs.String  = num2str(handles.phase_coeff{m,k},'%.2g');

% Set the edit field labels
handles.editMaxZ.String     = Zminmax;

if Action == 1
    handles.editWLmin.String    = minwl;
    handles.editWLmax.String    = maxwl;
    set(handles.plot_Ncontours,'String',Ncontours)
    set(handles.plot_colorrange,'String',plot_cRange)
    set(handles.plot_percentwhites,'String',plot_percentwhites)
end