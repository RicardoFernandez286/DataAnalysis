function dataStruct = ApplyChirpCorr(dataStruct,chirpSt,rootdir,type,k)

switch type
    case 'After fit'
        % Do nothing, the chirpSt is already given
        chirpFun    = chirpSt.chirpFun;
        Cfit        = chirpSt.Cfit;
    case 'Load fitted chirp'
        currdir = pwd;
        cd(rootdir);
        [fn,fp] = uigetfile('chirp*.mat','Select Chirp Correction File...','MultiSelect','off');
        cd(currdir)
        load([fp filesep fn]); %#ok<*LOAD>
end
%% Read from dataStruct
delays      = dataStruct.delays;
probeAxis   = dataStruct.cmprobe{k};

t0fit       = chirpFun(Cfit,probeAxis);

% Ndelays = length(t);
Npixels = length(probeAxis);

Zcorr   = dataStruct.corrdata{k};
Zraw    = dataStruct.rawsignal{k};
noise   = dataStruct.noise{k};

%% Interpolate data at each pixel
Zraw_new    = zeros(size(Zraw));
Zcorr_new   = zeros(size(Zcorr));
Noise_new   = zeros(size(noise));

for i=1:Npixels
    Zraw_new(:,i)   = interp1(delays-t0fit(i),Zraw(:,i),delays);
    Zcorr_new(:,i)  = interp1(delays-t0fit(i),Zcorr(:,i),delays);
    Noise_new(:,i)  = interp1(delays-t0fit(i),noise(:,i),delays);
end

Zraw_new    = Zraw_new(1:end-1,:);
Zcorr_new   = Zcorr_new(1:end-1,:);
Noise_new   = Noise_new(1:end-1,:);
delays      = delays(1:end-1);

switch dataStruct.rawcorr
    case 'CORRECTED'
        Z = dataStruct.corrdata;
    case 'RAW'
        Z = dataStruct.rawsignal;
end

%% Read the plot ranges
mintime     = min(delays);
maxtime     = max(delays);
minabs      = min(Z(:));
maxabs      = max(Z(:));
minwl       = min(probeAxis);
maxwl       = max(probeAxis);
Ncontours   = 40; % 40 contours by default is OK
plotranges  = [mintime maxtime minwl maxwl minabs maxabs Ncontours];

%% WRITE to dataStruct
% Write the main variables
dataStruct.delays       = delays;
dataStruct.cmprobe      = probeAxis;
dataStruct.rawsignal    = Zraw_new;
dataStruct.corrdata     = Zcorr_new;
dataStruct.noise        = Noise_new;
dataStruct.plotranges   = plotranges;
dataStruct.chirpCorr    = 1;

%% FOR DIAGNOSTICS
% %% Plot new data
% WHAT = Znew;
% maxPl = max(abs(WHAT(:)));
% minPl = -maxPl;
% ctrs  = linspace(minPl,maxPl,40);
% 
% fh = figure(1);
% clf(fh);
% fh.Color = 'w';
% 
% ax = axes('parent',fh);
% 
% contourf(ax,WL,t,WHAT,ctrs,'EdgeColor','flat');
% colormap(darkb2r(minPl,maxPl,40,2));
% caxis([minPl,maxPl]);
% colorbar;