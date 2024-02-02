function  dataStruct = load2DIRlab1(dataStruct,varargin)

% Description: This function loads all 2DIR data in the file format from the Lab 1 & Lab 4 MESS program.
% Usage: dataStruct = load2DIRlab1(dataStruct)
% Inputs:
%   dataStruct structure with fields:
%     datafilename
%     rootdir
%     
% Outputs:
%     cmprobe           (Double)
%     bins              (Double)
%     t2delays          (Double)
%     Ndelays           (Double)
%     Nspectra          (Double)
%     Ndatastates       (Double)
%     Nbins             (Double)
%     Nslowmod          (Double)
%     count             (Cell array)
%     probe             (Cell array)
%     reference         (Cell array)
%     interferogram     (Cell array)
%     signal			(Cell array)
%
% Ricardo Fernandez-Teran / 28.08.2019 / v5.0a

%% DEBUG
% datafilename = 'R3_mTiO2_Re1213_MeOH_4PM_2D_154022';
% rootdir = 'D:\Data\20180418';
% dataStruct.rootdir=rootdir;

%% HARDCODED Settings
autodetect_datatype = 1;
datatype = 'Raw'; % 'Raw' or 'Signal'

%% READ from dataStruct
if isempty(varargin)
    ShowWaitBar = true;
else
    ShowWaitBar = false;
end

datafilename    = dataStruct.datafilename;
rootdir         = dataStruct.rootdir;
Transient2D     = dataStruct.Transient2D;

%% Load all the necessary files after checking that they exist
filename        = [rootdir filesep datafilename filesep datafilename];

if exist([filename '_bins.csv'],'file') ~= 0
% Create progress bar and clear error string
if ShowWaitBar
    % dataStruct.WaitBar             = waitbar(0,'Loading data...');
    dataStruct.WBfigure                 = uifigure;
    dataStruct.WBfigure.Position(3:4)   = [405 175];
    dataStruct.WaitBar                  = uiprogressdlg(dataStruct.WBfigure,'Title','2D-IR data processing...','Message','Loading data...','Icon','info','ShowPercentage','on','Cancelable','on');
    dataStruct.ErrorText.String         = "";
end

% Load basic files
    cmprobe         = csvread([filename '_wavenumbers.csv']);
    bins            = csvread([filename '_bins.csv']);
    t2delays        = csvread([filename '_delays.csv']);
    
    % Determine whether the data is transient 2D or not, then read the transient delays
    if Transient2D == 1
        transient_delays    = csvread([filename '_transientDelays.csv']);
        t2delays            = [t2delays, transient_delays];
    end
        
    % TEMP dir story
    tempdir         = [rootdir filesep datafilename filesep 'temp'];
    tempdirOUT      = [rootdir filesep datafilename 'temp'];
    % Check if the temp dir is outside (Lab 1) or inside the data folder (Lab 4)
    if exist(tempdir,'dir') == 0 && exist(tempdirOUT,'dir') ~= 0
        tempdir     = tempdirOUT;
    end
    
% Continue loading other stuff
Nspectra        = csvread([filename '_Nspectra.csv']);
Ndatastates     = csvread([filename '_Ndatastates.csv']);
Nbins           = length(bins);
Ndelays         = length(t2delays);

% Interleaves and dummies
% Ninterleave     = length(csvread([filename '_interleaves.csv']));
Ninterleave     = 1;
Ndummies        = csvread([filename '_dummies.csv']);

% Read all the lines of the SlowModulation file
slowmodID       = fopen([filename '_slowModulation.csv']);
    % The SlowModulation file contains the No. of slowmods in the first line
    Nslowmod        = str2double(fgetl(slowmodID));
    % The rest of the lines contain the action matrices of all different modulations
    ModA            = strsplit(fgetl(slowmodID),':');
    ModB            = strsplit(fgetl(slowmodID),':');
    ModC            = strsplit(fgetl(slowmodID),':');
    ModSpec         = strsplit(fgetl(slowmodID),':');
    % In Lab 4, the piezo modulation is slow, so it is also included (in Lab 1 this would return -1)
    ModPiezo        = fgetl(slowmodID);
    % Close the file
fclose(slowmodID); clear slowmodID;

% Count the number of modulations for A,B,C and spectrometer:
N_modA          = sum(~isnan(str2num(ModA{2})),2);
N_modB          = sum(~isnan(str2num(ModB{2})),2);
N_modC          = sum(~isnan(str2num(ModC{2})),2);
N_modSpec       = sum(~isnan(str2num(ModSpec{2})),2);
% If any of them is zero (i.e. all NaN's), make it equal to 1 (i.e. 1 modulation = no modulation)
AllSlowMod                          = [N_modA N_modB N_modC N_modSpec];
AllSlowMod(AllSlowMod==0)           = 1;
N_modA          = AllSlowMod(1);
N_modB          = AllSlowMod(2);
N_modC          = AllSlowMod(3);
N_modSpec       = AllSlowMod(4);

% Get the spectrometer wavelengths (for spectrometer modulation only)
SpectrometerWL      = str2num(ModSpec{2});

% Count the number of 4-point modulation states (Lab 4 only)
if ModPiezo ~= -1 % if not at the end of the file
    ModPiezo    = strsplit(ModPiezo,':');
    N_4PM       = length(str2num(ModPiezo{2}))./(N_modA*N_modB*N_modC*N_modSpec);
else
    N_4PM = 1;
end

%% Load unfinished measurement from temp folder
if exist(tempdir,'dir') ~= 0
% Decide if we are previewing the data (i.e. the first scan didn't finish yet
    tempfilelist          = dir(tempdir);
    tempfilenames         = {tempfilelist.name}';
    % Select only "interferogram" files
    tempfilenames         = tempfilenames(contains(tempfilenames,"interferogram",'IgnoreCase',1));
    % Count how many delays have been finished
    finished_delays       = zeros(Ndelays,1);
    for i=1:Ndelays
        nfiles = length(tempfilenames(contains(tempfilenames,['de' num2str(i-1) '_'],'IgnoreCase',1)));
        if nfiles == Nslowmod*Ndatastates*Nspectra
            finished_delays(i)= 1;
        else
            finished_delays(i)= 0;
        end
    end
    % Count how many scans have been finished
    if length(tempfilenames) >= (Ndelays*Ndatastates*Nslowmod)
        Nscans = round(length(tempfilenames)/(Ndelays*Ndatastates*Nslowmod)); % Get the number of scans
    else
        Nscans = 0;
    end
    
    % Enter preview mode if no scans have been finished, otherwise continue normally
    if Nscans == 0
        preview_mode    = 1;
        Ndelays         = sum(finished_delays,1);
        t2delays        = t2delays(finished_delays == 1);
        autodetect_datatype = 0;
        datatype        = 'Raw';
        filename        = [tempdir filesep datafilename];
        dataStruct.ErrorText.String    = "Previewing data - Unfinished measurement!";
    else
        preview_mode    = 0; % Continue normally
        dataStruct.ErrorText.String    = "";
    end
    % If no delays have finished, warn the user
    if Ndelays == 0
        if ShowWaitBar
            warndlg('No delays completed yet!','No data');
            delete(dataStruct.WBfigure)
        end
        return
    end
else % Can't determine whether it's a preview or not - TEMP folder absent
    Ndelays         = length(t2delays);
    Nscans          = 1; % Set it to 1 so it's a valid number.
    preview_mode    = 0;
end

% Total number of states
totalstates     = Nspectra*Ndatastates*Nslowmod;
totalspectra    = totalstates*Ndelays;

%% Build a list of the files to be opened:
%   M= DELAYS, K=DATASTATES  - Will create a cell array of the form: t2 delays(rows) x states(columns)

% Reset the variables
endings={}; ShortEndings={}; n=0;

if preview_mode == 0 % IF AT LEAST ONE SCAN HAS FINISHED
    for q=0:Ndummies-1
        for p=0:Ninterleave-1
            for i=0:Nspectra-1
                for j=0:Nslowmod-1
                    for k=0:Ndatastates-1
                        n=n+1;
                        for m=0:Ndelays-1
                            endings{m+1,n}                          = ['_ds' num2str(k) '_sp' num2str(i) '_sm' num2str(j) '_de' num2str(m) '_in' num2str(p) '.csv'];
                            ShortEndings{m+1,ceil(n/Ndatastates)}   = ['_sp' num2str(i) '_sm' num2str(j) '_de' num2str(m) '_in' num2str(p) '.csv'];
                            endings_du{m+1,n}                       = ['_ds' num2str(k) '_sp' num2str(i) '_sm' num2str(j) '_de' num2str(m) '_in' num2str(p) '_du' num2str(q) '.csv'];
                            ShortEndings_du{m+1,ceil(n/Ndatastates)}= ['_sp' num2str(i) '_sm' num2str(j) '_de' num2str(m) '_in' num2str(p) '_du' num2str(q) '.csv'];
                        end
                    end
                end
            end
        end
    end

elseif preview_mode == 1 % IF THE MEASUREMENT IS STILL ONGOING AND THE FIRST SCAN DIDN'T FINISH YET
    for q=0:Ndummies-1
        for p=0:Ninterleave-1
            for i=0:Nspectra-1
                for j=0:Nslowmod-1
                    for k=0:Ndatastates-1
                        n=n+1;
                        for m=0:Ndelays-1
                            endings{m+1,n}      = ['_ds' num2str(k) '_sp' num2str(i) '_sm' num2str(j) '_de' num2str(m) '_in' num2str(p) '_0.csv'];
                            endings_du{m+1,n}   = ['_ds' num2str(k) '_sp' num2str(i) '_sm' num2str(j) '_de' num2str(m) '_in' num2str(p) '_du' num2str(q) '.csv'];
    %                       ShortEndings{m+1,round(n/Ndatastates)} = ['_sp' num2str(i) '_sm' num2str(j) '_de' num2str(m) '_in' num2str(p) '_0.csv'];
                        end
                    end
                end
            end
        end
    end
end
clear i j k m n;

% Verify that the RAW files exist (signal is a plus)
exist_count             = zeros(Ndelays,Ndatastates*Nslowmod*Ninterleave);
exist_probe             = zeros(Ndelays,Ndatastates*Nslowmod*Ninterleave);
exist_reference         = zeros(Ndelays,Ndatastates*Nslowmod*Ninterleave);
exist_interf            = zeros(Ndelays,Ndatastates*Nslowmod*Ninterleave);

% First, determine if new or old filename format
filelist        = dir([rootdir filesep datafilename]);
filenames       = {filelist.name};
filenames       = filenames(~contains(filenames,"dummies"));
find_dummies    = contains(filenames, "_du");
if sum(find_dummies(:)) > 0
    endings     = endings_du;
    ShortEndings= ShortEndings_du;
end

% Then, do the file existence test
for k=1:Ndatastates*Nslowmod*Ninterleave*Nspectra
    for m=1:Ndelays
      % For each dataset we need COUNTS, PROBE, REFERENCE & INTERFEROGRAM
        exist_count(m,k)            = exist([filename '_count' endings{m,k}],'file');
        exist_probe(m,k)            = exist([filename '_probe' endings{m,k}],'file');
        exist_reference(m,k)        = exist([filename '_reference' endings{m,k}],'file');
        exist_interf(m,k)           = exist([filename '_interferogram' endings{m,k}],'file');
    end
end

count_missing       = numel(exist_count(exist_count~=0))-numel(exist_count);
probe_missing       = numel(exist_probe(exist_probe~=0))-numel(exist_probe);
reference_missing   = numel(exist_reference(exist_reference~=0))-numel(exist_reference);
interf_missing      = numel(exist_interf(exist_interf~=0))-numel(exist_interf);

all_missing         = sum([count_missing;probe_missing;reference_missing;interf_missing]);

if all_missing == 0
% Figure out whether data is "RAW" or "SIGNAL" type (according to the MESS configuration when saving)
if autodetect_datatype == 1
    find_signal     = strfind(filenames, "signal");

    if sum([find_signal{:}],'omitnan') == 0
        datatype    = 'Raw';
    else
        datatype    = 'Signal';
    end
end

%% Calculate modulation type according to the number of datastates -NEEDS TO BE UPDATED
switch Ndatastates
    case 1
        Chopper = 'Chopper OFF';
    case 2
        Chopper = 'Chopper ON'; 
    case 4
        Chopper = 'Wobbler';
end

switch N_4PM
    case 1
        SlowMod = 'SlowMod OFF';
    case 4
        SlowMod = '4PM';
end

progress=0;
%% Open all the files according to the data type
switch datatype
case 'Raw'
    probe={}; reference={}; interferogram={}; count={};
    for k=1:Ndatastates*Nslowmod*Ninterleave*Nspectra
        for m=1:Ndelays
          % Update the Wait Bar
            if ShowWaitBar
                progress = progress+1;
                if dataStruct.WaitBar.CancelRequested
                    delete(dataStruct.WBfigure);
                    error('User aborted loading the data!');
                end
%                waitbar(progress/(Ndatastates*Nspectra*Nslowmod*Ndelays*Ninterleave*2),dataStruct.WaitBar,['Loading data (raw) - Reading file ' num2str(progress) ' of ' num2str(Ndatastates*Nslowmod*Ndelays*Ninterleave*Nspectra)]);
                dataStruct.WaitBar.Value    = progress/(Ndatastates*Nspectra*Nslowmod*Ndelays*Ninterleave*2);
                dataStruct.WaitBar.Message  = ['Loading data (raw) - Reading file ' num2str(progress) ' of ' num2str(Ndatastates*Nslowmod*Ndelays*Ninterleave*Nspectra)];    
            end
          % First load the counts
            count{m,k}          = csvread([filename '_count' endings{m,k}]);
            if cumsum(count{m,k}) == 0
                delete(dataStruct.WaitBar);
                error('2D-IR dataset incomplete - some delays have zero counts!')
            end
          % Then load probe, reference and interferogram
            probe{m,k}          = csvread([filename '_probe' endings{m,k}]);
            reference{m,k}      = csvread([filename '_reference' endings{m,k}]);
            interferogram{m,k}  = csvread([filename '_interferogram' endings{m,k}]);
          % Save only the first column of the interferogram
            interferogram{m,k}  = interferogram{m,k}(:,1);
          % Find the zero counts at the beginning & at the end
            start_NZ(m,k)       = find(count{m,k},1,'first');
            end_NZ(m,k)         = find(count{m,k},1,'last');
          % Remove the points with zero counts
            count{m,k}          = count{m,k}(start_NZ(m,k):end_NZ(m,k));
            probe{m,k}          = probe{m,k}(start_NZ(m,k):end_NZ(m,k),:);
            reference{m,k}      = reference{m,k}(start_NZ(m,k):end_NZ(m,k),:);
            interferogram{m,k}  = interferogram{m,k}(start_NZ(m,k):end_NZ(m,k));
            t1delays{m,k}       = bins(start_NZ(m,k):end_NZ(m,k),:);
          % Then, divide everything by the counts
            probe{m,k}          = probe{m,k}./count{m,k};
            reference{m,k}      = reference{m,k}./count{m,k};
            interferogram{m,k}  = interferogram{m,k}./count{m,k};
          % Next, divide by the average of each time trace to make the signal independent of the probe light spectrum
            probe{m,k}          = probe{m,k}./mean(probe{m,k},1);
            reference{m,k}      = reference{m,k}./mean(reference{m,k},1);
          % Now we can calculate the signal (in mOD)
            signal{m,k}         = -1000*log10(probe{m,k}./reference{m,k});
%             signal{m,k}         = -1000*(probe{m,k}./reference{m,k});
          % Store the size of the interferogram
            Int_size(m,k)       = length(interferogram{m,k});
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Finished loading everything

    % All data should have the same size
    for k=1:Ndatastates*Nslowmod*Ninterleave*Nspectra
        for m=1:Ndelays
            startcut            = Int_size(m,k)-min(Int_size);
            count{m,k}          = count{m,k}(startcut+1:end);
            probe{m,k}          = probe{m,k}(startcut+1:end,:);
            reference{m,k}      = reference{m,k}(startcut+1:end,:);
            interferogram{m,k}  = interferogram{m,k}(startcut+1:end);
            signal{m,k}         = signal{m,k}(startcut+1:end,:);
            t1delays{m,k}       = t1delays{m,k}(startcut+1:end,:);
        end
    end
    
%% Calculate the signal
    % Start with the interleaves - sum all interleaves
    if Ninterleave > 1
        for k=1:Ndatastates*Nslowmod*Nspectra
            for m=1:Ndelays
                for p=1:Ninterleave-1 
                    signal{m,k}         = signal{m,k}+signal{m,(k+p*Nslowmod*Ndatastates*Nspectra)};
                    interferogram{m,k}  = interferogram{m,k}+interferogram{m,(k+p*Nslowmod*Ndatastates*Nspectra)};
                end
                signal{m,k}             = signal{m,k}./Ninterleave;
                interferogram{m,k}      = interferogram{m,k}./Ninterleave;
            end
        end
        signal          = signal(:,1:Ndatastates*Nslowmod);
        interferogram   = interferogram(:,1:Ndatastates*Nslowmod);
    end
    
    % Continue with the slow modulation, it's either OFF or 4PM (polarisations and so on... TBI - not needed yet)
    switch SlowMod
    case 'SlowMod OFF' % This case applies for piezo random or off
        switch Chopper
        case 'Chopper ON' % Then we check the chopper
            % Calculate the signal by taking chopper ON - chopper OFF
            for m=1:Ndelays
                signal{m,1}=signal{m,1}-signal{m,2}; % Chopper ON - OFF
                interferogram{m,1}=interferogram{m,1}-interferogram{m,2};
            end
            % Reduce signal and interferogram to a single datastate (the calculated signal)
            signal=signal(:,1);
            interferogram=interferogram(:,1);
        case 'Chopper OFF'
            % If piezo is random/off and the chopper is off,
            % the signal is already there. Nothing to calculate
        case 'Wobbler'
            temp_sig=cell(Ndelays,Nslowmod*Nspectra); temp_int=cell(Ndelays,Nslowmod*Nspectra);
%             % Calculate the signal by taking the sum of all 4 states
%             for k=0:Nslowmod*Nspectra-1
%                 for m=1:Ndelays
%                   ds                = 1:Ndatastates;
% %                   temp_sig{m,k+1}   = sum(cat(3,signal{m,ds+k*Ndatastates}),3)/4;
%                   temp_sig{m,k+1}   = signal{m,1}-signal{m,2}+signal{m,3}-signal{m,4};
%                   temp_int{m,k+1}   = sum(cat(3,interferogram{m,ds+k*Ndatastates}),3)/4; 
%                 end
%             end
            
            % Calculate the signal by taking A-B+C-D
            for m=1:Ndelays
                signal{m,1}=signal{m,1}-signal{m,2}+signal{m,3}-signal{m,4}; % Chopper ON - OFF
                interferogram{m,1}=interferogram{m,1}-interferogram{m,2}+interferogram{m,3}-interferogram{m,4};
            end
%             % Save calculated signal and interferogram
%             signal                  = temp_sig;
%             interferogram           = temp_int;
        end
    case '4PM' % This case applies when we have four point modulation
        temp_sig=cell(Ndelays,Nslowmod*Nspectra); temp_int=cell(Ndelays,Nslowmod*Nspectra);
        % Calculate the signal for the four point modulation, but first check the chopper
        switch Chopper
        case 'Chopper ON'
            for m=1:Ndelays
                % Calculate the signal by taking chopper ON - chopper OFF for each SlowMod
                for k=0:Nslowmod-1
                    temp_sig{m,k+1}     = signal{m,2*(k+1)}-signal{m,2*k+1};
                    temp_int{m,k+1}     = interferogram{m,2*(k+1)}-interferogram{m,2*k+1};
                end
                % Then combine all datastates together for each delay
                temp_sig{m,1}           = sum(cat(3,temp_sig{m,:}),3)/4;
                temp_int{m,1}           = sum(cat(3,temp_int{m,:}),3)/4;
            end
            % Reduce signal and interferogram to a single datastate (the calculated signal)
            signal          = temp_sig(:,1);
            interferogram   = temp_int(:,1);
        case 'Chopper OFF'
            % Calculate the signal by averaging all four (modulation) points
            for m=1:Ndelays
                for k=1:Nslowmod
                    % Just combine all datastates together for each delay
                    temp_sig{m,1}       = sum(cat(3,signal{m,:}),3)/4;
                    temp_int{m,1}       = sum(cat(3,interferogram{m,:}),3)/4;
                end
            end
            signal          = temp_sig(:,1);
            interferogram   = temp_int(:,1);
        end
    end
    
case 'Signal'
    for k=1:Nslowmod*Ninterleave*Nspectra*Ndummies
    for m=1:Ndelays
        Ndatastates = 1;
        % Update the Wait Bar
        if ShowWaitBar
            progress = progress+1;
            if dataStruct.WaitBar.CancelRequested
                delete(dataStruct.WBfigure);
                error('User aborted loading the data!');
            end
    %         waitbar(progress/(Nslowmod*Ndelays*Ninterleave*Nspectra*Ndummies*2),dataStruct.WaitBar,['Loading data (signal) - Reading file ' num2str(progress) ' of ' num2str(Nslowmod*Ndelays*Ninterleave*Nspectra*Ndummies)]);
            dataStruct.WaitBar.Value    = progress/(Nslowmod*Ndelays*Ninterleave*Nspectra*Ndummies*2);
            dataStruct.WaitBar.Message  = ['Loading data (signal) - Reading file ' num2str(progress) ' of ' num2str(Nslowmod*Ndelays*Ninterleave*Nspectra*Ndummies)];
        end
      % Counts here are all 1's. Load them to ensure consistency...
        count{m,k}          = csvread([filename '_count' ShortEndings{m,k}]);
      % Load the interferogram and signal
        interferogram{m,k}  = csvread([filename '_interferogram' ShortEndings{m,k}]);
        signal{m,k}         = csvread([filename '_signal' ShortEndings{m,k}]);
      % Save only the first column of the interferogram
        interferogram{m,k}  = interferogram{m,k}(:,1);
      % Find the NaN in the arrays
        start_nNaN(m,k)     = find(~isnan(interferogram{m,k}),1,'first');
        end_nNaN(m,k)       = find(~isnan(interferogram{m,k}),1,'last');
      % Remove the points with zero counts (saved as NaNs by MESS)
        signal{m,k}         = signal{m,k}(start_nNaN(m,k):end_nNaN(m,k),:);
        interferogram{m,k}  = interferogram{m,k}(start_nNaN(m,k):end_nNaN(m,k));
        t1delays{m,k}       = bins(start_nNaN(m,k):end_nNaN(m,k),:);
      % Store the size of the interferogram
        Int_size(m,k)           = length(interferogram{m,k});
    end
    end
    for m=1:Ndelays
    % Cut the data to equal length
        startcut            = Int_size(m,k)-min(Int_size);
        count{m,k}          = count{m,k}(startcut+1:end);
        if strcmp(Chopper,'Wobbler')
            interferogram{m,k}  = interferogram{m,k}(startcut+1:end)/4;
            signal{m,k}         = signal{m,k}(startcut+1:end,:)/4;
        else
            interferogram{m,k}  = interferogram{m,k}(startcut+1:end);
            signal{m,k}         = signal{m,k}(startcut+1:end,:);
        end
        t1delays{m,k}       = t1delays{m,k}(startcut+1:end,:);
% %       %Normalize the interferogram
% %       interferogram{m,k}  = interferogram{m,k}./mean(interferogram{m,k});
    end
    
    % Continue with the slow modulation, it's either OFF or 4PM (polarisations and so on... TBI - not needed yet)
        switch SlowMod
        case '4PM' % This case applies when we have four point modulation
            temp_sig=cell(Ndelays,Ndummies);
            temp_int=cell(Ndelays,Ndummies);
            % Calculate the signal for the four point modulation, the chopper was already considered in MESS
            % Calculate the signal by averaging all four (modulation) points
                for k=0:Ndummies-1
                    for m=1:Ndelays
                        % Just combine all datastates together for each delay
                        sm              = 1:N_4PM;
                        temp_sig{m,k+1} = sum(cat(3,signal{m,sm+k*N_4PM}),3)/N_4PM;
                        temp_int{m,k+1} = sum(cat(3,interferogram{m,sm+k*N_4PM}),3)/N_4PM;
                    end
                end
                signal          = temp_sig;
                interferogram   = temp_int;
                Ndatastates     = 1; % Once we have the signal, set only one datastate
        end
end

%% Sort the data in increasing t2 (IR) delays
% These are needed
if size(t2delays,2) == 2
    [~,delay_index] = sort(t2delays(:,2),1);
else
    [~,delay_index] = sort(t2delays(:,1),1);
end
t2delays               = t2delays(delay_index,:);
signal                 = signal(delay_index,:);
interferogram          = interferogram(delay_index,:);

% These are also sorted for consistency
if strcmp(datatype,'Raw')
    count                  = count(delay_index,:);
    probe                  = probe(delay_index,:);
    reference              = reference(delay_index,:);
end


%% WRITE to dataStruct
dataStruct.isSimulation  = 0;
dataStruct.isShaper      = 0;
dataStruct.cmprobe       = cmprobe;
dataStruct.bins          = bins;
dataStruct.t2delays      = t2delays./1000; % in ps!
dataStruct.Ndelays       = Ndelays;
dataStruct.Nspectra      = Nspectra;
dataStruct.Ndummies      = Ndummies;
dataStruct.Ndatastates   = Ndatastates;
dataStruct.Nbins         = Nbins;
dataStruct.Nslowmod      = Nslowmod;
dataStruct.interferogram = interferogram;
dataStruct.signal		 = signal;
dataStruct.datatype      = datatype;
dataStruct.t1delays      = t1delays; % in fs!
dataStruct.Nscans        = Nscans;

% Clear 2DGC fit results
dataStruct.FitResults  = [];
dataStruct.t2_fitrange = [];
dataStruct.t2_fitdelays= [];
dataStruct.FitInput    = [];


if strcmp(datatype,'Raw')
    dataStruct.count     = count;
    dataStruct.probe     = probe;
    dataStruct.reference = reference;
end

else
    error('Not a valid 2D-IR dataset: some files missing!')
end
else
    error('Not a valid 2D-IR dataset: empty folder or corrupt data')
end