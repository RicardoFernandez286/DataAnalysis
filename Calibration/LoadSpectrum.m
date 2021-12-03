function data = LoadSpectrum(app,message)
% This subroutine loads a probe spectrum for further absorbance calculations
% Ricardo Fernandez-Teran / 09.11.2021 / v1.0a


[datafile, datapath, TYPE] = uigetfile({'*.2D','UoS TRIR';'*.dat','UniGE TA';'*.dat','UniGe nsTA'},message);
rootdir = app.rootdir;
if isempty(rootdir)
  rootdir = pwd;
end

if TYPE == 0
    return
end

currdir = pwd; 
cd(rootdir);

switch TYPE
    case 1 % Uni Sheffield TRIR
        DetSz   = [96 96 32 32]; % Sizes of [probe1 probe2 ref1 ref2] in pixels
        data    = cell(4,1);

        rawdata = readmatrix([datapath filesep datafile],'FileType','text','Delimiter','\t');
        for i=1:length(DetSz)
            idx =  (1:DetSz(i)) + sum(DetSz(1:i-1));
            data{i} = rawdata(idx);
        end      
        [~,fn,~]    = fileparts([datapath filesep datafile]);
        text        = readlines([datapath filesep fn '.LG']);
        g1_st       = strsplit(text{75},' '); % Det186 is 1st
        w1_st       = strsplit(text{77},' ');
        g2_st       = strsplit(text{67},' '); % Det185 is 2nd
        w2_st       = strsplit(text{69},' ');

        Gratings(1) = str2double(g1_st{end});
        CWL(1)      = str2double(w1_st{end});
        Gratings(2) = str2double(g2_st{end});
        CWL(2)      = str2double(w2_st{end});
        
        app.CAL_data.Gratings   = Gratings;
        app.CAL_data.CWL        = CWL;
        
    case 2 % UniGE TA
        rawdata = readmatrix([datapath filesep datafile],'FileType','text','Delimiter','\t');
        data{1} = rawdata(:,3);
                
    case 3 % UniGE nsTA
        rawdata = readmatrix([datapath filesep datafile],'FileType','text','Delimiter','\t');
        data{1} = rawdata(:,2);
        
end

app.CAL_data.CalType = TYPE;

cd(currdir);