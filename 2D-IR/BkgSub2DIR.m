function handles = BkgSub2DIR(handles)

%% READ from handles
phased_FFTZPsig     =   handles.phased_FFTZPsig;

%% Subtract the scattering background (if enabled)
for m=1:size()
    for k=1:size()
        if bkg_sub==1 && m~=1
            PROC_2D_DATA{m,k}       = phased_FFTZPsig{m,k}-phased_FFTZPsig{1,1};
        else
            PROC_2D_DATA{m,k}       = phased_FFTZPsig{m,k};
        end
    end
end