function data = LoadSpectrum(app,message,calType)
% This subroutine loads a probe spectrum for further absorbance calculations
% Ricardo Fernandez-Teran / 25.08.2023 / v2.0a

switch calType
    case 'UoS TRIR';                    TYPE = 1; uifilt = {'*.2D',calType};
    case 'UniGE TA';                    TYPE = 2; uifilt = {'*.dat',calType};
    case 'UniGE nsTA';                  TYPE = 3; uifilt = {'*.dat',calType};
    case 'UZH Lab 2';                   TYPE = 4; uifilt = {'*ds0_intensity*.csv',calType};
    case 'UniGE NIR-TA';                TYPE = 5; uifilt = {'*.dat',calType};
    case 'RAL LIFEtime (Absorbance)';   TYPE = 6; uifilt = {'*.csv',calType};
    case 'RAL LIFEtime (Intensity)';    TYPE = 7; uifilt = {'*cycle*.csv',calType};
    case 'UniGE TRIR (Intensity)';      TYPE = 8; uifilt = {'*ds0_intensity*.csv',calType};
    case 'UniGE TRIR (Absorbance)';     TYPE = 9; uifilt = {'*ds0_intensity*.csv',calType};
    otherwise;                          TYPE = 0;
end

[datafile, datapath] = uigetfile(uifilt,message);

rootdir = app.rootdir;
if isempty(rootdir)
  rootdir = pwd;
end

if TYPE == 0
    error('Invalid calibration data type.')
end

if datafile == 0
    error('No data selected!')
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
        
    case {2,5} % UniGE TA or NIR-TA
        rawdata = readmatrix([datapath filesep datafile],'FileType','text','Delimiter','\t');
        data{1} = rawdata(:,3);   
    
    case 3 % UniGE nsTA
        rawdata = readmatrix([datapath filesep datafile],'FileType','text','Delimiter','\t');
        data{1} = rawdata(:,2);
    
    case 4 % UZH Lab 2
        rawdata = readmatrix([datapath filesep datafile],'FileType','text','Delimiter','\t');
        data{1} = rawdata(:,1);
        
        [fp,~,~]    = fileparts(datapath);
        [fp,fn,~]   = fileparts(fp);
        text        = readlines([fp filesep fn filesep fn '_meta.txt']);
        g1_st       = strsplit(text{20},' ');
        w1_st       = strsplit(text{19},' ');
        
        Gratings(1) = str2double(g1_st{2});
        CWL(1)      = str2double(w1_st{3});
        
        app.CAL_data.Gratings   = Gratings;
        app.CAL_data.CWL        = CWL;
    case {6,7} % RAL LIFEtime (Absorbance)
        DetSz   = [128 128]; % Sizes of [probe1 probe2] in pixels -- no reference in LifeTime
        data    = cell(2,1);
        
        rawdata = readmatrix([datapath filesep datafile],'FileType','text','Delimiter',',');
        for i=1:length(DetSz)
            idx =  (1:DetSz(i)) + sum(DetSz(1:i-1));
            data{i} = rawdata(idx,2);
        end
        %         opts.Interpreter = 'tex';
        %         CWL = inputdlg({["Enter central wavenumbers for each detector";"(estimates are good, 0 to skip):";"Probe 1 [Left] (cm^{-1}):"];'Probe 2 [Right] (cm^{-1}):'},'Probe Calibration',[1,60],{'0','0'},opts);
        %         if isempty(CWL)
        %             error('Empty wavenumbers!')
        %         end
        %         CWL = str2double(CWL);
        %         if sum(CWL)==0
        %             error('Need at least one nonzero probe CWL!')
        %         end
        % 
        %         app.CAL_data.CWL = 1e7./CWL; % To NANOMETERS
    case {8,9}
        rawdata = readmatrix([datapath filesep datafile],'FileType','text','Delimiter','\t');
        data{1} = rawdata(:,1);
        
        [fp,~,~]    = fileparts(datapath);
        [fp,fn,~]   = fileparts(fp);
        text        = readlines([fp filesep fn filesep fn '_meta.txt']);
        
        g1_st       = strsplit(text{20},' ');
        w1_st       = strsplit(text{19},' ');
        
        if strcmp(w1_st{1},'Center wavelength') == 0
            g1_st       = strsplit(text{22},' ');
            w1_st       = strsplit(text{21},' ');
        end
        
        Gratings(1) = str2double(g1_st{2});
        CWL(1)      = str2double(w1_st{3});
        
        app.CAL_data.Gratings   = Gratings;
        app.CAL_data.CWL        = CWL;
end

app.CAL_data.CalType = TYPE;

cd(currdir);