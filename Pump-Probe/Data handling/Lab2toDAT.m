function Lab2toDAT(rootdir,datafilename)
% Usage: Lab2toDAT(rootdir,datafilename)
% Where both inputs are character arrays of the root directory and the 
% data name of the file to convert.
% 
% Will produce as an output a .dat file consisting of an array of the form:
%   (0)|  Wavenumbers
%    D |
%    e |
%    l |    DeltaAbs
%    a |    (in mOD)
%    y |
%    s |
% Ricardo Fernández-Terán
% v1.0 2017-07-24
    
if isemtpy(rootdir)
    rootdir=pwd;
end
cd(rootdir);

% Prepare the filenames for each file to read
delayfile = strcat(rootdir,filesep,datafilename,filesep,datafilename,'_delays.csv');
wavenumberfile = strcat(rootdir,filesep,datafilename,filesep,datafilename,'_wavenumbers.csv');
signalfile = strcat(rootdir,filesep,datafilename,filesep,datafilename,'_signal.csv');
noisefile = strcat(rootdir,filesep,datafilename,filesep,datafilename,'_signal_noise.csv');

% Read the files if the directory is correctly populated
delayfilech = char(delayfile);
wavenumberfilech = char(wavenumberfile);
signalfilech = char(signalfile);
noisefilech = char(noisefile);
existcondition = exist(delayfilech, 'file') && exist(wavenumberfilech, 'file') && exist(signalfilech, 'file') && exist(noisefilech,'file');
if exist(delayfilech, 'file') && exist(wavenumberfilech, 'file') && exist(signalfilech, 'file') && exist(noisefilech,'file')
    % Read the files
    delays = csvread(delayfile);
    cmprobe = csvread(wavenumberfile);
    rawsignal = csvread(signalfile);
    noise = csvread(noisefile);
    % Show by default the background subtraction limits
    % from the first data point till just before time zero
    k = 1;
    t = delays(k);
    while t < 0
        k = k+1;
        t = delays(k+1);
    end
    % Do the background subtraction and change status in handles.rawcorr
    bkg = mean(rawsignal(1:k,:));
    corrdata = rawsignal - bkg;
else
    assert(existcondition,'Selected directory is empty or files are corrupt!')
end
x=[0; delays];
y=transpose(cmprobe);
z=corrdata;
xyz=[x,[y;z]];
csvwrite([datafilename '.dat'],xyz);