% M= DELAYS, K=SCANS
%% Define critical variables
clear all
Ndelays = 2;
Ndatastates = 1;

% Initialise variables
    FFTinterferogram={}; binspecmax=[]; bininterfmax=[];
    bin=[]; difference={}; coeff={}; binzero=num2cell(-1*ones(Ndelays,Ndatastates));
    HeNe=2.11079; % Tau HeNe fs
    c_0=2.99792458e-5; % Speed of light in cm/fs
    
%% Read the inputs
    %Ndelays
    %Ndatastates
    %Nspectra
    %Nslowmod
    %Lastbin
    %BinZero 
    bkg_sub=0;
    filter_type='Median';
    filter_points=15;
    zeropad_factor=3;
    zeropad_next2k=1;
    zeropad_npoints=20;
    apodize_method='Cosine';
    apodize_gaussian=0;
    apodize_gaussiancoeff=20000;
    phase_method='Constant phase';
    phase_points=50;
    phase_factor=1;
    m=2;k=1;

    fmin=1850;d=6;
    fmax=fmin+31*d;
    cmprobe=fmin:d:fmax;
    freqmin=fmin;
    freqmax=fmax;
%% Read the data

% Read average interferograms (Pyro)
cd('D:\LABVIEW\2D ANALYSIS LAB4\Test data')
name=char("21.02.17.11'53'08 17.00");
interferogram{m,k}=dlmread([name ' _av pyro'],'\t');
coherence{m,k}=interferogram{m,k}(:,1);
interferogram{m,k}=interferogram{m,k}(:,2);

% Read average signal matrices (MCT)
signal{m,k}=dlmread([name ' _av int'],'\t');

% Take only the datapoints with a nonzero coherence delay
interferogram{m,k}=interferogram{m,k}(coherence{m,k} ~= 0);
signal{m,k}=signal{m,k}(coherence{m,k} ~= 0,:);

%% Preprocess the data
% Subtract slow background from interferogram and do FFT
switch filter_type
        case 'Median'
            interferogram{m,k}=interferogram{m,k}-mean(interferogram{m,k});
            interferogram{m,k}=interferogram{m,k}-medfilt1(interferogram{m,k},filter_points);
        case 'Mean'
            interferogram{m,k}=interferogram{m,k}-mean(interferogram{m,k});
end
FFTinterferogram{m,k}=fft(interferogram{m,k});

%% Find BinZero and phase before apodisation
% Find BinSpecMax and BinInterfMax
% Interferograms first, then use BinZero to phase the pixel data
    L=round(length(FFTinterferogram{m,k})/2);
    [~,binspecmax(m,k)]=max(abs(FFTinterferogram{m,k}(1:L)));
    [~,bininterfmax(m,k)]=max(abs(interferogram{m,k}(1:L)));
    % If binzero is not known (=-1) then try to guess it
    if binzero{m,k}==-1
        % The method will try shifting the interferogram until the point where the phase is flat
        for p=1:100
            bin(p)=p+binspecmax(m,k);
            tempFFT{m,k}=fft(flip(circshift(interferogram{m,k},-bin(p))));
            % Unwrap the phases to remove discontinuities greater than ±?
            tempPhase{m,k}=unwrap(angle(tempFFT{m,k}));
            difference{m,k}(p)=(tempPhase{m,k}(binspecmax(m,k)+10))-(tempPhase{m,k}(binspecmax(m,k)-10));
        end
        % Then, it will fit a straight line that tells how much the slope is changing
        coeff{m,k}=polyfit(bin,difference{m,k},1);
        % Afterwards, it will find the zero crossing point from the fitted coefficients
        binzero{m,k}=round(-(coeff{m,k}(2)/coeff{m,k}(1)));
    end

% Get the interferogram phase
    % Define variables
    PhasedInt={}; FFTPhasedInt={}; FFTintFlatPhase={}; phase={};
    % Shift the interferogram to phase it correctly (negative shift!)
    PhasedInt{m,k}=circshift(interferogram{m,k},-binzero{m,k});
    FFTPhasedInt{m,k}=fft(PhasedInt{m,k});
    % Unwrap the phases to remove discontinuities greater than ±?
    FFTintFlatPhase{m,k}=unwrap(angle(FFTPhasedInt{m,k}));
    % Get the phase (in radians)
    phase{m,k}=mod(FFTintFlatPhase{m,k}(binspecmax(m,k)),pi);
% Cut the interferogram and pixel data (remove BinZero points from the end) / 1:end-binzero{m,k}-1
    interferogram{m,k}=interferogram{m,k}(1:end-binzero{m,k});
    signal{m,k}=signal{m,k}(1:end-binzero{m,k},:);
    Nbins_cut{m,k}=length(interferogram{m,k});
    
%% Zeropad, apodize and phase the data
% Get the points for zeropadding
if apodize_method=="Box"
    mean_interf{m,k}=mean(interferogram{m,k}(end-zeropad_npoints:end),1);
    mean_sig{m,k}=mean(signal{m,k}(end-zeropad_npoints:end,:),1);
else
    mean_interf{m,k}=0;
    mean_sig{m,k}=zeros(1,length(cmprobe));
end

% Apodize the data by the selected method
    % Cosine
    switch apodize_method
        case 'Box'
            cos_exp=0;
        case 'Cosine'
            cos_exp=1;
        case 'Cosine squared'
            cos_exp=2;
        case 'Cosine cube'
            cos_exp=3;
    end
    cosine{m,k} = zeros(Nbins_cut{m,k},1);
    box{m,k} = zeros(Nbins_cut{m,k},1);
    for q=1:Nbins_cut{m,k}
        cosine{m,k}(q)=(cos(pi*(q-binzero{m,k})/(2*(Nbins_cut{m,k}-binzero{m,k})))^cos_exp);
        box{m,k}(q)=(heaviside(q-binzero{m,k})-heaviside(q-Nbins_cut{m,k}));
    end
    cosine{m,k}(end)=0;

    % Gaussian
    gaussian{m,k} = zeros(Nbins_cut{m,k},1);
    apodize_function{m,k} = zeros(Nbins_cut{m,k},1);
    if apodize_gaussian==1
        for q=1:Nbins_cut{m,k}
            gaussian{m,k}(q)=exp(-((q-binzero{m,k}))^2/(apodize_gaussiancoeff*Nbins_cut{m,k}));
        end
        apodize_function{m,k}=cosine{m,k}.*gaussian{m,k};
    else
        apodize_function{m,k}=cosine{m,k};
    end

% Multiply the interferogram by the apodization function
    apo_interferogram{m,k}=interferogram{m,k}.*apodize_function{m,k};
% Multiply the signal by the apodization function
    apo_signal{m,k}=signal{m,k}.*apodize_function{m,k};

% Build the matrices
    N_newpoints{m,k}=zeropad_factor*Nbins_cut{m,k};
    if zeropad_next2k==1
        N_newpoints{m,k}=2^nextpow2(N_newpoints{m,k});
    end
    N_newpoints{m,k}=N_newpoints{m,k}-Nbins_cut{m,k};
    pad_interf{m,k}=ones(N_newpoints{m,k},1)*mean_interf{m,k};
    pad_signal{m,k}=ones(N_newpoints{m,k},length(cmprobe)).*mean_sig{m,k};

% Do the zeropadding
    zeropad_interf{m,k}=[apo_interferogram{m,k};pad_interf{m,k}];
    zeropad_signal{m,k}=[apo_signal{m,k};pad_signal{m,k}];

% Find the phase with the known BinZero
    [~,binspecmax(m,k)]=max(abs(fft(zeropad_interf{m,k})));
    [~,bininterfmax(m,k)]=max(abs(zeropad_interf{m,k}));

% Shift the interferogram to BinZero and get the phase
    phased_ZPint{m,k}=circshift(zeropad_interf{m,k},-binzero{m,k});
    phased_ZPsig{m,k}=circshift(zeropad_signal{m,k},-binzero{m,k});
    Npoints=length(phased_ZPint{m,k});

% Calculate the FFT and phase
    FFT_ZPint{m,k}=fft(phased_ZPint{m,k},[],1);
    FFT_ZPsig{m,k}=fft(phased_ZPsig{m,k},[],1);
    ZP_phase{m,k}=mod(unwrap(angle(FFT_ZPint{m,k})),pi);
    
% Fit the phase to the whole interferogram
    points=transpose([floor(binspecmax(m,k)-phase_points/2):1:floor(binspecmax(m,k)+phase_points/2)]);
    switch phase_method
        case 'Constant phase'
            fittedPhase{m,k}=ones(length(ZP_phase{m,k}),1)*mean(ZP_phase{m,k}(points));
        case 'Linear phase'
            phase_coeff{m,k}=polyfit(points,ZP_phase{m,k}(points),1);
            fittedPhase{m,k}=transpose(polyval(phase_coeff{m,k},1:1:Npoints));
        case 'Quadratic phase'
            phase_coeff{m,k}=polyfit(points,ZP_phase{m,k}(points),2);
            fittedPhase{m,k}=transpose(polyval(phase_coeff{m,k},1:1:Npoints));
    end

% Phase the 2D data
    phase{m,k}=fittedPhase{m,k}(binspecmax(m,k));
    phasingterm{m,k}=exp(-1i*fittedPhase{m,k}*phase_factor);
    phased_FFTZPsig{m,k}=real(FFT_ZPsig{m,k}.*phasingterm{m,k});
    phased_FFTZPint{m,k}=real(FFT_ZPint{m,k}.*phasingterm{m,k});
    magnitude_FFTsig{m,k}=abs(FFT_ZPsig{m,k});
    
% Calculate the pump frequency axis from interferogram by using HeNe counting
    freq_axis=transpose([1:1:Npoints])/(Npoints*HeNe*c_0);

%% PROBE AXIS CALIBRATION
% Get a list of the maxima of the interferogram (first half)
    P=floor(Npoints/200);
    range=(binspecmax(m,k)-P):(binspecmax(m,k)+P);
    [~,maxindex]=max(magnitude_FFTsig{m,k}(range,:));
    plot(freq_axis(range(maxindex)),'-o')
    Npixels = size(magnitude_FFTsig{m,k},2);
    pixels=transpose(1:1:Npixels);
    
    freq_coeff{m,k}=polyfit(pixels,freq_axis(range(maxindex)),3);
    freq_fit=transpose(polyval(freq_coeff{m,k},pixels));
    plot(freq_fit)
% Set the limits from the selected frequency ranges (default min/max of cmprobe)
    [bin_limits,~]= findClosestId2Val(freq_axis,[freqmin,freqmax]);
    sel_bins=bin_limits(1):1:bin_limits(2);
    sel_freq=freq_axis(sel_bins);

%% BACKGROUND SUBTRACTION 
% Subtract the background (if enabled)
    if bkg_sub==1 && m~=1
        PROC_2D_DATA{m,k}=phased_FFTZPsig{m,k}(sel_bins,:)-phased_FFTZPsig{1,1}(sel_bins,:);
    else
        PROC_2D_DATA{m,k}=phased_FFTZPsig{m,k}(sel_bins,:);
    end

    
% plot_2DIR
% plot(freq_axis,phased_FFTZPint{m,k})
