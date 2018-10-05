% function LoadTRES()
%% Get the data folder and filenames
% Ask for the folder
folder      = uigetdir();
cd(folder);
% Get list of files
files       = dir(folder);
% Remove directories
files       = files(~[files.isdir]);
% Remove TRES_data.dat (if it's there)
files(strcmp({files.name},'TRES_data.dat')) = [];
% Remove IRF.txt (if it's there)
files(strcmp({files.name},'IRF.txt'))       = [];

%% Read the files
xy  = dlmread(files(1).name,'\t',10,0);
z   =[0 sscanf(files(1).name,strcat('Decay [','%f','nm]'))];
WB  = waitbar(1/length(files),['Loading data - File 1 of ' num2str(length(files)-1)]);
for i = 2:length(files)
    % Waitbar
    waitbar(i/(length(files)),WB,['Loading data - File ' num2str(i) ' of ' num2str(length(files)-1)]);
    % Read the file
    currentfile = files(i).name;
    data = dlmread(currentfile,'\t',10,0);
    xy = [xy data(:,2)];
    % In this case, Z is the wavelength (nm)
    currentz = sscanf(currentfile,strcat('Decay [','%f','nm]'));
    z = [z currentz];
end
delete(WB)

%% Process the data
%%%%%%%%% Time Calibration = ns/ch
timecal = 0.02743484;
xy(:,1) = xy(:,1)*timecal;
%%%%%%%%% Time shift (find global intensity maximum and set as zero).
[~,MaxWL] = max(max(xy(1:end,2:end)));
[~,MaxT] = max(xy(1:end,MaxWL));
tshift = xy(MaxT,1);
xy(:,1) = xy(:,1) - tshift;
%%%%%%%%% Concatenate matrices (WL / data)
xyz = vertcat(z,xy);

%% Plot the data (ask user)
% Ask
choice = questdlg('Plot the data?','Plot','2D plot','3D plot','No, thank you','No, thank you');

% Plot
switch choice
    case '2D plot'
        fh = figure();
        contourf(xyz(2:end,1),xyz(1,2:end),transpose(xyz(2:end,2:end)),20,'LineStyle','-','LineColor','flat');
    case '3D plot'
        fh = figure();
        s=surf(xyz(2:end,1),xyz(1,2:end),transpose(xyz(2:end,2:end))); xlim([-2,10]);
        s.EdgeColor='none';
        s.FaceColor='interp';
end
        fh.Color = 'White';
        colorbar;
        title('Time-resolved emission spectra');
        xlabel('Time (ns)','FontSize',12,'FontWeight','bold');
        ylabel('Wavelength (nm)','FontSize',12,'FontWeight','bold');

%% SAVE the data
savedata = questdlg('Save the data?','Save?','Yes','No','Yes');
switch savedata
    case 'Yes'
        filename = [pwd filesep 'TRES_data.dat'];
        dlmwrite(filename,xyz);
        helpdlg(['Data saved in ' filename],'Data saved!')
end

%% Do SVD
yy = xy(:,2:end);
[svd_time,svals,svd_WL] = svd(yy,'econ');

svd_WL      = transpose(svd_WL*svals);
svd_time    = transpose(svd_time);
svals       = diag(svals);
