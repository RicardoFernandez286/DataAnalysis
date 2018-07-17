clear all

%% READ inputs
% datedir     = 'D:\Data\20180314';
% samplename  = 'H2O_FeCN_dil_225817';
% blankname   = 'H2O_230247';

datedir     = 'D:\Data\20180511';
samplename  = 'R06_mTiO2_Re1213_MeOH_FTIR_MCT_174832';
blankname   = 'blank_FTIR_141549';

datatype_blank      = 't0 fast'; % '2DIR' or 't0 fast'
datatype_sample     = 't0 fast';

BL_corr_type='All'; % 'None', 'Moving', 'All', 'Offset'
threshold   = 20; % Spectral threshold in % of the maximum
filt_order  = 4;

% For 2DIR, which delay to use?
delay_spl   = 1;
delay_blk   = 1;

Ndatastates = 2;
numpyro     = 1;
useblank    = 1;

filter_type = 'Mean';
filter_points= 20;

box         = 0;
cos_exp     = 5;
gaussian_apo= 0;
triangle_apo= 0;

gaussian_exp= 5000000;
phase_points= 20;

do_zeropad  = 1;
phase_method= 'Constant';

PlotAll     = 0;

% If we are not using the blank, then just copy the sample name
if useblank == 0
    blankname = samplename;
end

% Physical constants
HeNe            = 2.11079;          % HeNe period (fs)
c_0             = 2.99792458e-5;    % Speed of light in cm/fs

%% Read the files
blk_count={}; blk_int={};
spl_count={}; spl_int={};
blank={};
sample = {};

for m=1:Ndatastates
    % Read the blank
    switch datatype_blank
        case '2DIR'
            blk_count{m}= csvread([datedir filesep blankname filesep blankname '_count_ds' num2str(m-1) '_sp0_sm0_de' num2str(delay_blk-1) '_in0.csv']);
            blk_int{m}  = csvread([datedir filesep blankname filesep blankname '_interferogram_ds' num2str(m-1) '_sp0_sm0_de' num2str(delay_blk-1) '_in0.csv']);
        case 't0 fast'
            blk_count{m}= csvread([datedir filesep blankname filesep blankname '_sp0_ds' num2str(m-1) '_count_0.csv']);
            blk_int{m}  = csvread([datedir filesep blankname filesep blankname '_sp0_ds' num2str(m-1) '_interferogram_0.csv'])';
    end
    
    % Read the sample
    switch datatype_sample
        case '2DIR'
            spl_count{m}= csvread([datedir filesep samplename filesep samplename '_count_ds' num2str(m-1) '_sp0_sm0_de' num2str(delay_spl-1) '_in0.csv']);
            spl_int{m}  = csvread([datedir filesep samplename filesep samplename '_interferogram_ds' num2str(m-1) '_sp0_sm0_de' num2str(delay_spl-1) '_in0.csv']);
        case 't0 fast'
            spl_count{m}= csvread([datedir filesep samplename filesep samplename '_sp0_ds' num2str(m-1) '_count_0.csv']);
            spl_int{m}  = csvread([datedir filesep samplename filesep samplename '_sp0_ds' num2str(m-1) '_interferogram_0.csv'])';
    end
    % Match the sizes of the arrays (i.e. 2 interferograms - 2 count columns (the same))
    if size(blk_int{m},2) > 1
        blk_count{m}=blk_count{m}*ones(1,size(blk_int{m},2));
    end
    if size(spl_int{m},2) > 1
        spl_count{m}=spl_count{m}*ones(1,size(spl_int{m},2));
    end
    
    % Find the zero counts at the beginning & at the end (SAMPLE)
        start_NZ_spl(m)     = find(spl_count{m}(:,1),1,'first');
        end_NZ_spl(m)       = find(spl_count{m}(:,1),1,'last');
    % Remove the points with zero counts (BLANK)
        spl_count{m}        = spl_count{m}(start_NZ_spl(m):end_NZ_spl(m),:);
        spl_int{m}          = spl_int{m}(start_NZ_spl(m):end_NZ_spl(m),:);

    % Find the zero counts at the beginning & at the end (BLANK)
        start_NZ_blk(m)     = find(blk_count{m}(:,1),1,'first');
        end_NZ_blk(m)       = find(blk_count{m}(:,1),1,'last');
    % Remove the points with zero counts (BLANK)
        blk_count{m}        = blk_count{m}(start_NZ_blk(m):end_NZ_blk(m),:);
        blk_int{m}          = blk_int{m}(start_NZ_blk(m):end_NZ_blk(m),:);
        
    % Divide everything by the counts to get actual volts
    blank{m}    = blk_int{m} ./ blk_count{m};
    sample{m}   = spl_int{m} ./ spl_count{m};

    for i=1:2
        switch filter_type
            case 'Median'
                sample{m}(:,i)      = -(sample{m}(:,i)-mean(sample{m}(:,i)));
                sample{m}(:,i)      = sample{m}(:,i)-medfilt1(sample{m}(:,i),filter_points);
                blank{m}(:,i)       = -(blank{m}(:,i)-mean(blank{m}(:,i)));
                blank{m}(:,i)       = blank{m}(:,i)-medfilt1(blank{m}(:,i),filter_points);
            case 'Mean'
                sample{m}(:,i)      = -(sample{m}(:,i) - mean(sample{m}(:,i)));
                blank{m}(:,i)       = -(blank{m}(:,i) - mean(blank{m}(:,i)));
        end
     end
end

%% Process the data
% Calculate interferograms (CHOPPER ON - CHOPPER OFF)
sample_int      = sample{1} - sample{2};
blank_int       = blank{1} - blank{2};

FFT_spl         = fft(sample_int);
FFT_blk         = fft(blank_int);

Npoints_spl     = length(sample_int);
Npoints_blk     = length(blank_int);

% Find BinZero in order to apodise the interferogram, then zeropad it 
    %%%% Find BinSpecMax for the SAMPLE [do it using the INTERFEROMETER pyro -> (~,1)]
    [~,bininterfmax_spl]    = max(sample_int(:,1));
    [~,binspecmax_spl]      = max(abs(FFT_spl(20:floor(Npoints_spl/2),1)));
    binspecmax_spl          = binspecmax_spl + 20;
    % The method will try shifting the interferogram until the point where the phase is flat
        for p=1:100
            bin(p)              = p+binspecmax_spl;
            tempFFT             = fft(flip(circshift(sample_int(:,1),-bin(p))));
            % Unwrap the phases to remove discontinuities greater than +-pi
            tempPhase           = unwrap(angle(tempFFT));
            difference(p)       = (tempPhase(binspecmax_spl+10))-(tempPhase(binspecmax_spl-10));
        end
        % Then, it will fit a straight line that tells how much the slope is changing
        coeff           = polyfit(bin,difference,1);
        % Afterwards, it will find the zero crossing point from the fitted coefficients
        binzero_spl     = round(-(coeff(2)/coeff(1)));

    %%%% Find BinSpecMax and BinInterfMax for the BLANK [do it using the INTERFEROMETER pyro -> (~,1)]
    [~,bininterfmax_blk]    = max(blank_int(:,1));
    [~,binspecmax_blk]      = max(abs(FFT_blk(20:floor(Npoints_blk/2),1)));
    binspecmax_blk          = binspecmax_blk + 20;
    % The method will try shifting the interferogram until the point where the phase is flat
        for p=1:100
            bin(p)              = p+binspecmax_blk;
            tempFFT             = fft(flip(circshift(blank_int(:,1),-bin(p))));
            % Unwrap the phases to remove discontinuities greater than +-pi
            tempPhase           = unwrap(angle(tempFFT));
            difference(p)       = (tempPhase(binspecmax_blk+10))-(tempPhase(binspecmax_blk-10));
        end
        % Then, it will fit a straight line that tells how much the slope is changing
        coeff           = polyfit(bin,difference,1);
        % Afterwards, it will find the zero crossing point from the fitted coefficients
        binzero_blk     = round(-(coeff(2)/coeff(1)));

%% Apodise the data
    %%% Apodise the SAMPLE
    cosine_spl                  = zeros(Npoints_spl,1);
    gaussian_spl                = zeros(Npoints_spl,1);
    box_spl                     = zeros(Npoints_spl,1);
    for q=1:Npoints_spl
        cosine_spl(q)           = cos(pi*(q-binzero_spl)/(2*(Npoints_spl-binzero_spl)))^cos_exp;
        gaussian_spl(q)         = exp(-((q-binzero_spl)/(gaussian_exp./Npoints_spl))^2);
        box_spl(q)              = heaviside(q-binzero_spl);
    end
    if gaussian_apo == 1
        apo_sample              = sample_int.*cosine_spl.*gaussian_spl;
    elseif gaussian_apo == 0
        apo_sample              = sample_int.*cosine_spl;
    end
    
    %%% Apodise the BLANK
    cosine_blk                  = zeros(Npoints_blk,1);
    gaussian_blk                = zeros(Npoints_blk,1);
    box_blk                     = zeros(Npoints_blk,1);
    for q=1:Npoints_blk
        cosine_blk(q)           = cos(pi*(q-binzero_blk)/(2*(Npoints_blk-binzero_blk)))^cos_exp;
        gaussian_blk(q)         = exp(-((q-binzero_blk)/(gaussian_exp./Npoints_blk))^2);
        box_blk(q)              = heaviside(q-binzero_blk);
    end
    if gaussian_apo == 1
        apo_blank               = blank_int.*cosine_blk.*gaussian_blk;
    else
        apo_blank               = blank_int.*cosine_blk;
    end
    
    %%% Apply box (throw away negative times)
    if box == 1
        apo_sample              = apo_sample.*box_spl;
        apo_blank               = apo_blank.*box_blk;
    end
    
%% Zeropad the data by a factor 2, then to the next 2^k    
    NpointsFFT      = 2^nextpow2(2*max(Npoints_blk,Npoints_spl));  
    N_newpoints_spl = NpointsFFT-Npoints_spl;
    N_newpoints_blk = NpointsFFT-Npoints_blk;
    pad_matrix_spl  = zeros(N_newpoints_spl,2);
    pad_matrix_blk  = zeros(N_newpoints_blk,2);
    
    if do_zeropad == 1
        ZPsample_int= [apo_sample;pad_matrix_spl];
        ZPblank_int = [apo_blank;pad_matrix_blk];
    else
        ZPsample_int= apo_sample;
        ZPblank_int = apo_blank;
        NpointsFFT  = length(apo_sample);
    end

% Calculate the frequency axis
    Resolution  = 1/(NpointsFFT*HeNe*c_0);
    freq        = ((1:1:NpointsFFT)-1)'.*Resolution;
       
%% Phase the data
    % Shift to start at BinZero
    shifted_ZPspl       = circshift(ZPsample_int,-binzero_spl);
    shifted_ZPblk       = circshift(ZPblank_int,-binzero_blk);
    % Do FFT
    FFT_ZPsample        = fft(shifted_ZPspl);
    FFT_ZPblank         = fft(shifted_ZPblk);

    
    % Calculate the points for the phase fit
    shift_fitphase              = round(phase_points/Resolution);
    points_spl                  = transpose((binspecmax_spl-shift_fitphase):1:(binspecmax_spl+shift_fitphase));
    points_blk                  = transpose((binspecmax_blk-shift_fitphase):1:(binspecmax_blk+shift_fitphase));
    
    for i=1:2
        % Unwrap the phases
        ZP_phase_spl        = angle(FFT_ZPsample(:,i));
        ZP_phase_blk        = angle(FFT_ZPblank(:,i));
        % Do the phase fit
        switch phase_method
            case 'Constant'
                fittedPhase_spl     = ones(length(ZP_phase_spl),1)*mean(ZP_phase_spl(points_spl));
                fittedPhase_blk     = ones(length(ZP_phase_blk),1)*mean(ZP_phase_blk(points_blk));
            case 'Linear'
                phase_coeff_spl     = polyfit(points_spl,ZP_phase_spl(points_spl),1);
                fittedPhase_spl     = transpose(polyval(phase_coeff_spl,1:1:NpointsFFT));
                phase_coeff_blk     = polyfit(points_blk,ZP_phase_blk(points_blk),1);
                fittedPhase_blk     = transpose(polyval(phase_coeff_blk,1:1:NpointsFFT));
            case 'Quadratic'
                phase_coeff_spl     = polyfit(points_spl,ZP_phase_spl(points_spl),2);
                fittedPhase_spl     = transpose(polyval(phase_coeff_spl,1:1:NpointsFFT));
                phase_coeff_blk     = polyfit(points_blk,ZP_phase_blk(points_blk),2);
                fittedPhase_blk     = transpose(polyval(phase_coeff_blk,1:1:NpointsFFT));
            case 'Cubic'
                phase_coeff_spl     = polyfit(points_spl,ZP_phase_spl(points_spl),3);
                fittedPhase_spl     = transpose(polyval(phase_coeff_spl,1:1:NpointsFFT));
                phase_coeff_blk     = polyfit(points_blk,ZP_phase_blk(points_blk),3);
                fittedPhase_blk     = transpose(polyval(phase_coeff_blk,1:1:NpointsFFT));
            case 'No fit'
                fittedPhase_spl     = ZP_phase_spl;
                fittedPhase_blk     = ZP_phase_blk;
        end
        % Phase the FT data
        Sample_spectrum(:,i)        = real(FFT_ZPsample(:,i).*exp(-1i*fittedPhase_spl));
        Blank_spectrum(:,i)         = real(FFT_ZPblank(:,i).*exp(-1i*fittedPhase_blk));
    end
%% Calculate absorption spectra
if numpyro == 2
    sample_spec = -log10(Sample_spectrum(:,2)./Sample_spectrum(:,1));
    blank_spec  = -log10(Blank_spectrum(:,2)./Blank_spectrum(:,1));
elseif numpyro == 1
    sample_spec = -log10(Sample_spectrum(:,2));
    blank_spec  = -log10(Blank_spectrum(:,2));
end

if useblank == 1
    Absorbance = real(sample_spec - blank_spec);
else
    Absorbance = real(sample_spec);
end

% Select useful sepctral range
accepted_points = (abs(Sample_spectrum(1:round(NpointsFFT/2),1)) >= (threshold/100).*max(abs(Sample_spectrum(1:round(NpointsFFT/2),1))));
spec_start      = find(accepted_points,1,'first');
spec_end        = find(accepted_points,1,'last');

Absorbance      = Absorbance(spec_start:spec_end);
freq            = freq(spec_start:spec_end);

% Baseline correction 'None', 'Moving', 'All', 'Offset'
switch BL_corr_type
    case 'All'
        Corr_Abs = Absorbance - medfilt1(Absorbance,filt_order,[],'truncate');
        Corr_Abs = 1000*(Corr_Abs - min(Corr_Abs));
    case 'Moving'
        Corr_Abs = 1000.*(Absorbance - medfilt1(Absorbance,filt_order,[],'truncate'));
    case 'Offset'
        Corr_Abs = 1000*(Absorbance - min(Absorbance));
    case 'None'
        Corr_Abs = 1000.*Absorbance;
end

%% Plot
subplot(2,1,1)
plot(freq,Corr_Abs,'-b','DisplayName','Subtracted spectrum','LineWidth',2);

if PlotAll == 1
    hold on
    plot(freq,1000.*sample_spec(spec_start:spec_end),'DisplayName','Sample spectrum');
    plot(freq,1000.*blank_spec(spec_start:spec_end),'DisplayName','Blank spectrum');
    hold off
    legend
end
rl = refline(0,0);
rl.Color = [0.5 0.5 0.5];
rl.Annotation.LegendInformation.IconDisplayStyle = 'off';

xlim([freq(1) freq(end)]);
ylim('auto')
title('In-situ FTIR spectrum','FontWeight','bold');
ylabel('Absorbance (mOD)','FontWeight','bold');
xlabel('Wavenumbers (cm^{-1})','FontWeight','bold');

subplot(2,1,2)
plot(ZPsample_int)