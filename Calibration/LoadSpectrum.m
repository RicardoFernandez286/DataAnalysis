function data = LoadSpectrum(app,message)
DetSz   = [96 96 32 32]; % Sizes of [probe1 probe2 ref1 ref2] in pixels
data    = cell(4,1);
rootdir = app.rootdir;
if isempty(rootdir)
    rootdir = pwd;
end
currdir = pwd;

cd(rootdir);

[datafile, datapath] = uigetfile('*.2D',message);

rawdata = readmatrix([datapath filesep datafile],'FileType','text','Delimiter','\t');

for i=1:length(DetSz)
    idx =  [1:DetSz(i)] + sum(DetSz(1:i-1));
    data{i} = rawdata(idx);
end

cd(currdir);
end