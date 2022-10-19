function SplitCalibration
% This script will merge two calibration files
% (i.e. take the first detector and the second detector from different
% files and write it back into a single file).

saveULTRA = 1;

[f1,p1] = uigetfile('*.csv','Select DETECTOR 1 (LHS) calibration','MultiSelect','off');
[f2,p2] = uigetfile('*.csv','Select DETECTOR 2 (RHS) calibration','MultiSelect','off');

% Read files
cm1 = readmatrix([p1 f1]);
cm2 = readmatrix([p2 f2]);

% Check the lenghts are equal (!)
if length(cm1) ~= length(cm2)
    error('Unequal calibrated probe lenghts! Check')
end

%% Merge and write
% Take the 1st half from f1 and the 2nd half from f2
Npix = round(length(cm1)/2);

cm_merged = [cm1(1:Npix); cm2(Npix+(1:Npix))];

pW = uigetdir(pwd,'Save MERGED calibration file to...');

if pW == 0
    warning('No directory selected, aborting!')
    return
end

try 
    writematrix(cm_merged,[pW 'CalibratedProbe.csv']);
    if saveULTRA == 1
        writematrix(cm_merged,[(1:(2*Npix))' [pW sprintf('CalibProbe_%s.csv',datestr(now,'yyyymmdd-HHMMSS'))]]);
    end
    warning('Merged calibration file saved successfully!')
catch err
    warning('Error writing merged calibration files')
end
