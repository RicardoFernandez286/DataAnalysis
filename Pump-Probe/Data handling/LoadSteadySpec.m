function dataStruct = LoadSteadySpec(dataStruct,varargin)
if ~isempty(varargin)
    SettingsPath    = varargin{1};
else
    SettingsPath    = 'C:\GUIoptions.txt';
end

% Get settings file
fileexist = exist(SettingsPath,'file');
if fileexist == 2
    try 
        defaultSSdir = readParam('defaultSSdir',SettingsPath);
    catch
        warning('Default SS directory not defined in GUIoptions.txt');
        defaultSSdir = pwd;
    end
    if exist(defaultSSdir,'dir') == 0
        defaultSSdir = pwd;
    end
    old_dir = pwd;
    cd(defaultSSdir);
    [FileName,PathName,FilterIndex] = uigetfile(...
    {'*.dat;*.csv;*.dpt','ASCII spectra (*.dat,*.csv,*.dpt)';'*.??','Bruker OPUS files (*.#)'}, ...
    'Select the FTIR/UV-Vis spectrum to load...');
else
    [FileName,PathName,FilterIndex] = uigetfile(...
    {'*.dat;*.csv;*.dpt','ASCII spectra (*.dat,*.csv,*.dpt)';'*.??','Bruker OPUS files (*.#)'}, ...
    'Select the FTIR/UV-Vis spectrum to load...');
    old_dir = pwd;
end

cd(old_dir);

if FileName ~= 0
    SteadySpecfile = [PathName FileName];
    % Load it according to format: Either OPUS or dat/csv/pdt
    % X axis in Wavenumbers assumed
    switch FilterIndex
        case 2
            try
                [y,x,~] = ImportOpus(SteadySpecfile,'RatioAbsorptionChanged');
            catch
                try
                    [y,x,~] = ImportOpus(SteadySpecfile);
                catch
                    return
                end
            end
        case 1
            try
                SSdata  = readmatrix(SteadySpecfile);
            catch err
                try
                    SSdata  = readmatrix(SteadySpecfile,'Delimiter',',');
                catch err
                    return
                end
            end
            x       = SSdata(:,1);
            y       = SSdata(:,2);
    end
    % Convert Y data to mOD
    y = 1000.*y;
    % Store data in dataStruct
    dataStruct.SSx          = x;
    dataStruct.SSy          = y;
    dataStruct.SSloaded   = 1;
else
    dataStruct.SSloaded   = 0;
end