function batch2DGC_fit(rootfolder,fitType,varargin)
%% Define startup variables
% rootfolder  = 'D:\Ricardo Data\switchdrive\Ph.D. UZH\MANUSCRIPTS\9) 2D IR distance - na\Latest Simulations\Dimer_distance1';
% fitType     = 'Test'; % 'Ricardo' or 'Andrea'

%% Build list of (sub)folders, which contain the data
folderslist = dir(rootfolder);
foldernames = {folderslist.name}';
foldernames = foldernames([folderslist.isdir]);
foldernames = foldernames(3:end);
Ndatafiles  = length(foldernames);

% Read more options from varargin:
%   varargin{1} = No. of t2 delays to skip
%   varargin{2} = Path to a text file containing a list of the datasets to fit
if ~isempty(varargin)
    Nskip   = varargin{1};
    disp([datestr(now,-1) ': Fitting every ' num2str(Nskip) ' points in t2.']);
    
    if length(varargin) > 1
        disp([datestr(now,-1) ': Reading directory list...']);
        fid = fopen(varargin{2});
        j=0;
        while ~feof(fid)
            j=j+1;
            imp_datapath{j} = fgetl(fid);
        end
        imp_datapath = imp_datapath';
        Ndatafiles = length(imp_datapath);
        fclose(fid);
        disp([datestr(now,-1) ': Found ' num2str(Ndatafiles) ' datasets in list.']);
    end
else
    Nskip   = [];
end

if Ndatafiles == 0
    disp([datestr(now,-1) ': Nothing to do here! Bye :)']);
    exit
end

%% Start parallel pool (if none exists)
% if isempty(gcp('nocreate'))
%    parpool;
% end

%% Loop through the folders: Load, process and fit all the spectra
for i=1:Ndatafiles
    if length(varargin) > 1
%         [folderpath,foldernames{i}] = fileparts(imp_datapath{i});
        imp_parts   = strsplit(imp_datapath{i},filesep);
        folderpath  = [];
        for s=1:length(imp_parts)-1
            folderpath = [folderpath filesep imp_parts{s}];
        end
        foldernames{i} = imp_parts{end};
        rootfolder  = ['/home/group' filesep folderpath];
    end

    %%% Determine the data type
    if exist([rootfolder filesep foldernames{i} filesep 'spec_2D.dat'],'file') ~= 0
        dataType                = 'Simulation';
    elseif exist([rootfolder filesep foldernames{i} filesep foldernames{i} '_bins.csv'],'file') ~= 0
        dataType                = 'Experiment';
    else
        disp([datestr(now,-1) ': ' 'Folder ' foldernames{i} ' does not contain a valid 2D-IR dataset']);
        continue % Not a valid dataset, try with the the next one
    end
    
    %%% Load the data
    dataStruct.rootdir      = rootfolder;
    dataStruct.datafilename = foldernames{i};
    load([fileparts(mfilename('fullpath')) filesep 'defaultGUI.mat']);
    
    switch dataType
        case 'Simulation'
            disp([datestr(now,-1) ': ' 'Loading dataset ' foldernames{i} ' [Simulation]']);
            dataStruct = load2DIRsimu(dataStruct,'NoWaitBar');
        case 'Experiment'
            disp([datestr(now,-1) ': ' 'Loading dataset ' foldernames{i} ' [Experiment]']);
            dataStruct.Transient2D = 0;
            dataStruct = load2DIRlab1(dataStruct,'NoWaitBar');
    end
    
    %%% Get the (default) fit parameters
    switch fitType
        case 'Ricardo'
%             cut_data = 'Use probe axis';
            cut_data = 'VCH2-HF';
			fitparameters =    ...
				[{'2210'}
				 {'0'   }
				 {'80'  }
				 {'12'  }
				 {'8'   }
				 {'1'   }
				 {'Diag'}];
% 				[{'1998'}    {'2098'}    {'1998' }    {'2098' }    {'1967'}    {'2018'}
% 				 {'0'   }    {'0'   }    {'2098' }    {'1998' }    {'0'   }    {'0'   }
% 				 {'25'  }    {'24'  }    {'26'   }    {'24'   }    {'20'  }    {'31'  }
% 				 {'12'  }    {'17'  }    {'12'   }    {'17'   }    {'12'  }    {'12'  }
% 				 {'8'   }    {'15'  }    {'15'   }    {'8'    }    {'9'   }    {'9'   }
% 				 {'1'   }    {'1'   }    {'1'    }    {'1'    }    {'1'   }    {'1'   }
% 				 {'Diag'}    {'Diag'}    {'Xpeak'}    {'Xpeak'}    {'Diag'}    {'Diag'}];
            %t2_fitrange = [0.25 max(dataStruct.t2delays)];
			t2_fitrange = [0.25 40];
            equal_SxSy  = 0;
            diffSyfor12 = 1;
            writepath   = '/home/ricfer/FitResults/Solvation/VCH2-HF';
        otherwise
            cut_data    = 'Use probe axis';
            fitparameters =    ...      
           [{'1979'}    {'2028'}    {'1979' }    {'2028'  }
            {'0'   }    {'0'   }    {'2028' }    {'1979'  }
            {'10'  }    {'10'  }    {'-1'   }    {'-1'    }
            {'8'   }    {'8'   }    {'-1'   }    {'-1'    }
            {'8'   }    {'8'   }    {'-1'   }    {'-1'    }
            {'1'   }    {'1'   }    {'0h'   }    {'0h'    }
            {'Diag'}    {'Diag'}    {'Xpeak'}    {'Xpeak'}];
            t2_fitrange = [9, 11];
            equal_SxSy  = 1;
            diffSyfor12 = 0;
            writepath   = '/home/ricfer/FitResults';
    end
    
    %%% Check we haven't fitted it before (tiple-check twice...)
%     if exist([rootfolder filesep foldernames{i} '_FIT_RESULTS.mat'],'file') ~= 0
%         disp([datestr(now,-1) ': ' 'Dataset ' foldernames{i} ' has been fitted before! Delete/rename the old fit file before trying to fit again.']);
%         continue
%     end
    
    if exist([writepath filesep foldernames{i} '_FIT_RESULTS.mat'],'file') ~= 0
        disp([datestr(now,-1) ': ' 'Dataset ' foldernames{i} ' has been fitted before! Delete/rename the old fit file before trying to fit again.']);
        continue
    end

    %%% Process the data (experiment only)
    switch dataType
        case 'Experiment'
            disp([datestr(now,-1) ': ' 'Processing ' foldernames{i} ' ...']);
            app.I2D_AutocalibrateprobeaxisCheckBox.Value = 0;
            dataStruct = process2DIR(app,dataStruct,0,'NoWaitBar');
    end

    %%% Do the fit
    try
        disp([datestr(now,-1) ': ' 'Starting fit of ' foldernames{i} ' ...']);
        [dataStruct,~] = Gaussian2D_analysis(app,dataStruct,cut_data,fitparameters,t2_fitrange,equal_SxSy,diffSyfor12,writepath,Nskip);
    catch
        disp([datestr(now,-1) ': ' 'Error fitting ' foldernames{i} ' ...']);
        continue
    end
end

disp([datestr(now,-1) ': ' 'Everything done!']);
% exit