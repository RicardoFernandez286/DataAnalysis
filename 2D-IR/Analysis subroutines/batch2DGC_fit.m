function batch2DGC_fit(rootfolder,fitType,varargin)
%% Define startup variables
rootfolder  = 'D:\Ricardo Data\switchdrive\Ph.D. UZH\MANUSCRIPTS\9) 2D IR distance - na\Latest Simulations\Dimer_distance1';
fitType     = 'Test'; % 'Ricardo' or 'Andrea'

%% Start parallel pool (if none exists)
if isempty(gcp('nocreate'))
   parpool;
end

%% Build list of (sub)folders, which contain the data
folderlist  = dir(rootfolder);
foldernames = {folderlist.name}';
foldernames = foldernames([folderlist.isdir]);
foldernames = foldernames(3:end);
Ndatafiles  = length(foldernames);

if ~isempty(varargin)
    Nskip   = varargin{1};
else
    Nskip   = [];
end

%% Loop through the folders: Load, process and fit all the spectra
for i=1:Ndatafiles
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
            cut_data = 'Re1213 VET';
            fitparameters =    ...
           [{'1979'}    {'2028'}    {'1979' }    {'2028'  }
            {'10'  }    {'10'  }    {'-2028'}    {'-1979' }
            {'8'   }    {'8'   }    {'-1'   }    {'-1'    }
            {'8'   }    {'8'   }    {'-1'   }    {'-1'    }
            {'1'   }    {'1'   }    {'0h'   }    {'0h'    }
            {'Diag'}    {'Diag'}    {'Xpeak'}    {'Xpeak'}];
            t2_fitrange = [0.1 max(dataStruct.t2delays)];
            equal_SxSy  = 0;
            diffSyfor12 = 1;
            writepath   = '/home/ricfer/FitResults';
        case 'Andrea'
            cut_data = 'Use probe axis';
            fitparameters = ...
           [{'1980'}    {'2060'}    {'1980' }    {'2060'  }
            {'25'  }    {'25'  }    {'-2060'}    {'-1980' }
            {'15'  }    {'15'  }    {'-1'   }    {'-1'    }
            {'10'  }    {'10'  }    {'-1'   }    {'-1'    }
            {'1'   }    {'1'   }    {'0h'   }    {'0h'    }
            {'Diag'}    {'Diag'}    {'Xpeak'}    {'Xpeak'}];
            t2_fitrange = [0.5 max(dataStruct.t2delays)];
            equal_SxSy  = 1;
            diffSyfor12 = 0;
            writepath   = '/home/apasti/FitResults';
        otherwise
            cut_data    = 'Use probe axis';
            fitparameters =    ...
           [{'1979'}    {'2028'}    {'1979' }    {'2028'  }
            {'11'  }    {'11'  }    {'-2028'}    {'-1979' }
            {'10'  }    {'10'  }    {'-1'   }    {'-1'    }
            {'10'  }    {'10'  }    {'-1'   }    {'-1'    }
            {'1'   }    {'1'   }    {'0h'   }    {'0h'    }
            {'Diag'}    {'Diag'}    {'Xpeak'}    {'Xpeak'}];
            t2_fitrange = [9, 11];
            equal_SxSy  = 1;
            diffSyfor12 = 0;
            writepath   = '/home/ricfer/FitResults';
    end
    
    %%% Check we haven't fitted it before (tiple-check twice...)
    if exist([rootfolder filesep foldernames{i} '_FIT_RESULTS.mat'],'file') ~= 0
        disp([datestr(now,-1) ': ' 'Dataset ' foldernames{i} ' has been fitted before! Delete/rename the old fit file before trying to fit again.']);
        continue
    end
    
    if exist([writepath filesep foldernames{i} '_FIT_RESULTS.mat'],'file') ~= 0
        disp([datestr(now,-1) ': ' 'Dataset ' foldernames{i} ' has been fitted before! Delete/rename the old fit file before trying to fit again.']);
        continue
    end

    %%% Process the data (experiment only)
    switch dataType
        case 'Experiment'
            disp([datestr(now,-1) ': ' 'Processing ' foldernames{i} ' ...']);
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

% exit