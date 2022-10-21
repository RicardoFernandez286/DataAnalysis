function SplitCalibration
% This script will merge two calibration files
% (i.e. take the first detector and the second detector from different
% files and write it back into a single file).

saveULTRA = 1;

[f1,p1] = uigetfile('*.csv','Select Calibration File');

% Read files
cm_merged = readmatrix([p1 f1]);

% Check the lenghts are equal (!)
if mod(length(cm_merged),2) ~= 0
    error('Odd number of pixels!')
end

%% Merge and write
% Take the 1st half from f1 and the 2nd half from f2
Npix = round(length(cm_merged)/2);

pW = uigetdir(pwd,'Save SPLIT calibration files to...');

if pW == 0
    warning('No directory selected, aborting!')
    return
end

try 
    writematrix([(1:(Npix))' cm_merged(1:Npix)], [pW filesep sprintf('CalibProbe_%s_LHS.csv',datestr(now,'yyyymmdd-HHMMSS'))]);
    writematrix([(Npix + (1:(Npix)))' cm_merged(Npix+(1:Npix))], [pW filesep sprintf('CalibProbe_%s_RHS.csv',datestr(now,'yyyymmdd-HHMMSS'))]);
    disp('Split calibration files saved successfully!')
catch err
    warning('Error writing merged calibration files')
end
