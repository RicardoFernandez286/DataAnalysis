function dataStruct = LoadSteadyIR(dataStruct,varargin)
if ~isempty(varargin)
    SettingsPath    = varargin{1};
else
    SettingsPath    = 'C:\GUIoptions.txt';
end

% Get settings file
fileexist = exist(SettingsPath,'file');
if fileexist == 2
    try 
        defaultIRdir = readParam('defaultIRdir',SettingsPath);
    catch
        warning('Default IR directory not defined in GUIoptions.txt');
        defaultIRdir = pwd;
    end
    if exist(defaultIRdir,'dir') == 0
        defaultIRdir = pwd;
    end
    old_dir = pwd;
    [FileName,PathName,FilterIndex] = uigetfile(...
    {'*.dat;*.csv;*.dpt','ASCII spectra (*.dat,*.csv,*.dpt)';'*.??','Bruker OPUS files (*.#)'}, ...
    'Select the FTIR spectrum to load...');
else
    [FileName,PathName,FilterIndex] = uigetfile(...
    {'*.dat;*.csv;*.dpt','ASCII spectra (*.dat,*.csv,*.dpt)';'*.??','Bruker OPUS files (*.#)'}, ...
    'Select the FTIR spectrum to load...');
    old_dir = pwd;
end

cd(old_dir);

if FileName ~= 0
    SteadyIRfile = [PathName FileName];
    % Load it according to format: Either OPUS or dat/csv/pdt
    % X axis in Wavenumbers assumed
    switch FilterIndex
        case 2
            try
                [y,x,~] = ImportOpus(SteadyIRfile,'RatioAbsorptionChanged');
            catch
                try
                    [y,x,~] = ImportOpus(SteadyIRfile);
                catch
                    return
                end
            end
        case 1
            try
                IRdata  = readmatrix(SteadyIRfile);
            catch err
                try
                    IRdata  = readmatrix(SteadyIRfile,'Delimiter',',');
                catch err
                    return
                end
            end
            x       = IRdata(:,1);
            y       = IRdata(:,2);
    end
    % Convert Y data to mOD
    y = 1000.*y;
    % Store data in dataStruct
    dataStruct.IRx          = x;
    dataStruct.IRy          = y;
    dataStruct.FTIRloaded   = 1;
else
    dataStruct.FTIRloaded   = 0;
end