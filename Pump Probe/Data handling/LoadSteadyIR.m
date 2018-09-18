function handles = LoadSteadyIR(handles)
% Get file
fileexist = exist('C:\GUIoptions.txt','file');
if fileexist == 2
    try 
        defaultIRdir = readParam('defaultIRdir','C:\GUIoptions.txt');
    catch
        warning('Default IR directory not defined in GUIoptions.txt');
        defaultIRdir = pwd;
    end
    cd(defaultIRdir);
    [FileName,PathName,FilterIndex] = uigetfile(...
    {'*.dat;*.csv;*.dpt','ASCII spectra (*.dat,*.csv,*.dpt)';'*.??','Bruker OPUS files (*.#)'}, ...
    'Select the FTIR spectrum to load...');
else
    [FileName,PathName,FilterIndex] = uigetfile(...
    {'*.dat;*.csv;*.dpt','ASCII spectra (*.dat,*.csv,*.dpt)';'*.??','Bruker OPUS files (*.#)'}, ...
    'Select the FTIR spectrum to load...');
end
if FileName ~= 0
    SteadyIRfile = [PathName FileName];
    % Load it according to format: Either OPUS or dat/csv/pdt
    % X axis in Wavenumbers assumed
    switch FilterIndex
        case 2
            try
                [y,x,~] = ImportOpus(SteadyIRfile,'RatioAbsorption');
            catch
                try
                    [y,x,~] = ImportOpus(SteadyIRfile);
                catch
                    return
                end
            end
        case 1
            IRdata  =  csvread(SteadyIRfile);
            x       = IRdata(:,1);
            y       = IRdata(:,2);
    end
    % Convert Y data to mOD
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