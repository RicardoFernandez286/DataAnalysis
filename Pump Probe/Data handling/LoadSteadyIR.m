function handles = LoadSteadyIR(handles)
% Get file
fileexist = exist('C:\GUIoptions.txt','file');
if fileexist == 2
    defaultIRdir = readParam('defaultIRdir','C:\GUIoptions.txt');
    cd(defaultIRdir);
    [FileName,PathName,FilterIndex] = uigetfile(...
    {'*.dat;*.csv;*.dpt','ASCII spectra (*.dat,*.csv,*.dpt)';'*.??','Bruker OPUS files (*.??)'}, ...
    'Select the FTIR spectrum to load...');
else
    [FileName,PathName,FilterIndex] = uigetfile(...
    {'*.dat;*.csv;*.dpt','ASCII spectra (*.dat,*.csv,*.dpt)';'*.??','Bruker OPUS files (*.??)'}, ...
    'Select the FTIR spectrum to load...');
end
if FileName ~= 0
    SteadyIRfile = [PathName FileName];
    % Load it according to format: Either OPUS or dat/csv/pdt
    % X axis in Wavenumbers assumed
    switch FilterIndex
        case 2
            [y,x,~] = ImportOpus(SteadyIRfile,'RatioAbsorption');
        case 1
            IRdata =  csvread(SteadyIRfile);
            x = IRdata(:,1);
            y = IRdata(:,2);
    end
    % Y data in mODS
    y = 1000.*y;
    % Store data in Handles
    handles.IRx = x;
    handles.IRy = y;
    handles.FTIRloaded.Value = 1;
    handles.ShowFTIR.Visible = 'On';
    handles.EraseFTIR.Visible = 'On';
else
    handles.IRloaded = 0;
end