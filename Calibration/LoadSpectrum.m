function data = LoadSpectrum(app,message)

[datafile, datapath, idx] = uigetfile({'*.2D','UoS TRIR';'*.dat','UniGE TA'},message);
rootdir = app.rootdir;
if isempty(rootdir)
  rootdir = pwd;
end
currdir = pwd; 
cd(rootdir);

switch idx
    case 1 % Uni Sheffield TRIR
        DetSz   = [96 96 32 32]; % Sizes of [probe1 probe2 ref1 ref2] in pixels
        data    = cell(4,1);

        rawdata = readmatrix([datapath filesep datafile],'FileType','text','Delimiter','\t');
        for i=1:length(DetSz)
            idx =  (1:DetSz(i)) + sum(DetSz(1:i-1));
            data{i} = rawdata(idx);
        end
        app.CAL_data.CalType = 1;
    
    case 2 % UniGE TA
        rawdata = readmatrix([datapath filesep datafile],'FileType','text','Delimiter','\t');
        data{1} = rawdata(:,3);
        app.CAL_data.CalType = 2;
end

cd(currdir);