function CheckFitResults(plotDelay)
%% Define startup variables
% rootfolder  = 'D:\Ricardo Data\switchdrive\Ph.D. UZH\MANUSCRIPTS\9) 2D IR distance - na\Latest Simulations\Dimer_distance1';
% fitfolder   = [rootfolder filesep 'FitResults_new'];
rootfolder  = 'D:\Ricardo Data\switchdrive\Ph.D. UZH\RESULTS\2D-IR\Lab 4\Second round\20181002';
fitfolder   = rootfolder;

plotDelay       = 20;
multiPlot       = 0;
plotResiduals   = 0;

%% Build list of (sub)folders, which contain the data
folderlist  = dir(rootfolder);
foldernames = {folderlist.name}';
foldernames = foldernames([folderlist.isdir]);
foldernames = foldernames(3:end);
Ndatafiles  = length(foldernames);

[idx,tf]    = listdlg('ListString',foldernames);

if tf == 0
    return
end

Nplots = length(idx);

if multiPlot == 1
    fh = figure(1);
    fh.Color = [1 1 1];
    fh.Units = 'normalized';
    fh.OuterPosition = [0 0 1 1];
    ax1 = axes(fh);
    n_sub = ceil(sqrt(Nplots));
    m_sub = ceil(sqrt(n_sub^2-Nplots));
end

for k=1:Nplots
    i=idx(k);
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
            load([fileparts(mfilename('fullpath')) filesep 'defaultGUI.mat']);
            disp([datestr(now,-1) ': ' 'Processing ' foldernames{i} ' ...']);
            app.I2D_AutocalibrateprobeaxisCheckBox.Value = 1;
            dataStruct = process2DIR(app,dataStruct,0,'NoWaitBar');
    end

    %%% Load the fit results .MAT file
    if exist([fitfolder filesep foldernames{i} '_FIT_RESULTS.mat'],'file') ~= 0
        if multiPlot == 1
            ax = subplot(n_sub,n_sub,k);
        else
            % Create figure    
            fh = figure(k);
            ax = axes(fh);
            ax.Visible='Off';
            fh.Color = [1 1 1];

            % Make uniform, consistent format
            ax.FontSize = 12;
            ax.LineWidth = 1;
            ax.TickLength = [0.015 0.035];
            ax.Visible='On';
        end
        
        % Read default plot options
        load([fileparts(mfilename('fullpath')) filesep 'defaultGUI.mat']);
        load([fileparts(mfilename('fullpath')) filesep 'defaultPlotOptions.mat']);
        plotOptions.popdelay        = findClosestId2Val(dataStruct.t2delays,plotDelay);
        plotOptions.plotFitResults  = 1;
        plotOptions.plotResidualsFit= plotResiduals;
        dataStruct.Transient2D      = 0;
        % Read the fit results
        load([fitfolder filesep foldernames{i} '_FIT_RESULTS.mat']);
        dataStruct.FitResults   = FitResults;
        dataStruct.t2_fitrange  = minmax(dataStruct.t2delays');
        dataStruct.t2_fitdelays = t2delays;
        if exist('input','var') == 1
            input_st = load([fitfolder filesep foldernames{i} '_FIT_RESULTS.mat'],'input');
            dataStruct.FitInput     = input_st.input;
        else    
            dataStruct.FitInput     = input_st;
        end
        % Do the plot
        ContourPlot_2DIR(plotOptions,dataStruct,ax);
        title(ax,foldernames{i},'interpreter','none');
    else
        disp([datestr(now,-1) ': ' 'No fit file found for dataset ' foldernames{i} '(!)']);
        continue
    end

end
