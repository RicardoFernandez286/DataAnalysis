function Photoisomerisation()
% Calculate and plot photoisomerisation difference spectra from .txt files
% CURRENT VERSION v1.0

% Get the directory
ConvertDir = uigetdir();
% If the user selected a directory
if ConvertDir ~= 0
% Ask the user what to do
    choice = questdlg('Convert all files from OPUS to CSV?', ...
        'Convert files?', ...
        'Yes','No','No');
    switch choice
        case 'No'
            convert = 0;
        case 'Yes'
            convert = 1;
    end
    % Select only files that are NOT *.dpt, *.dat or *.csv
    list = dir(ConvertDir); list(1:2,:) = [];
    allfiles = transpose({list.name});
    I = contains(allfiles,["dpt","csv","dat"]);
    files = allfiles(I==0);
    cd(ConvertDir);
    switch convert
       case 1
        prompt = {'Enter the time between succesive spectra:','Enter time units:'};
        dlg_title = 'Time between spectra';
        num_lines = 1;
        defaultans = {'2', 's'};
        timespec = inputdlg(prompt,dlg_title,num_lines,defaultans);
        time = str2num(timespec{1});
        units = timespec{2};
        % Read all files and write everything in a matrix
        xvector = []; header = [];
            data=dlmread(char(files(1)));
            x=data(:,1);
            ymatrix = zeros(length(x),length(files));
            ymatrix(:,1) = data(:,2);
        for i=2:length(files)
            currfile = char(files(i));
            data=dlmread(currfile);
            ymatrix(:,i) = data(:,2);
        end
        header = time .* (0:1:length(files)-1);
        ymatrix = vertcat(header,ymatrix);
        xvector = vertcat(0,x);
        csvwrite('ALLDATA.dat',[xvector,ymatrix])
        uiwait(msgbox('Conversion completed!', 'Success!','warn'));
        % Ask the user what to do next
        Struct.Interpreter = 'tex';
        Struct.WindowStyle = 'modal';
        Struct.Default = 'No';
        delta = questdlg('Calculate \DeltaAbs matrix?', ...
            'Calculate difference spectra?', ...
            'Yes, please','Yes and plot','No',Struct);
        i=0; j=0;
        switch delta
            case 'Yes, please'
                calc=1; makeplot=0;
            case 'Yes and plot'
                calc=1; makeplot=1;
        end
        if calc == 1
                % Calculate Delta_abs matrix
                prompt = ['Indicate the number of spectra to average for the baseline. ' 'Remember, you only have ' num2str(size(ymatrix,2)) ' spectra, so choose wisely'];
                n = inputdlg(prompt,'Baseline averaging',1,{'5'});
                n = str2num(n{:});
                deltaymatrix = []; header=[];
                % Average the first N spectra for the baseline
                baseline = mean(ymatrix(:,1:n),2);
                for i=1:size(ymatrix,2)
                    deltaymatrix = [deltaymatrix,ymatrix(:,i)-baseline];
                end
                header = time .* (1-n:1:length(files)-n);
                deltaymatrix = vertcat(header,deltaymatrix(2:end,:));
                % Save it
                csvwrite('DELTA_ABS.dat',[xvector,deltaymatrix])
                % Ask for a reduced range and save it
                WL = inputdlg('Define the spectral window to save:','Define spectral window',1,{'210 700'});
                WL = str2num(WL{:}); minWL=WL(1); maxWL=WL(2);
                name = ['DELTA_ABS_' num2str(minWL) ' to ' num2str(maxWL) 'nm.dat'];
                deltaymatrixselect=[];
                [minI, ~] = findClosestId2Val(xvector,minWL);
                [maxI, ~] = findClosestId2Val(xvector,maxWL);
                deltaymatrixselect = vertcat(header,deltaymatrix(minI:maxI,:));
                csvwrite(name,[vertcat(0,xvector(minI:maxI)),deltaymatrixselect])
            if makeplot == 1
                % Plot the stuff with a rainbow and nice axes
                figure; hold on;
                cmap = colormap(jet(size(deltaymatrix,2)));
                for i=1:size(deltaymatrix,2)
                    plot(xvector(2:end),deltaymatrix(2:end,i),'Color',cmap(i,:));
                end
                hline = refline(0,0); hline.Color = [0.5 0.5 0.5];
                set(gca,'FontSize',14)
                xlabel('Wavelength (nm)','FontSize',14,'FontWeight','bold')
                ylabel('\DeltaAbs (OD)','FontSize',14,'FontWeight','bold')
                xlim([xvector(2),xvector(end)]); caxis([0 deltaymatrix(1,end)]); hcb=colorbar;
                n = round(deltaymatrix(1,end)./6);
                hcb.YTick = 0:n:deltaymatrix(1,end);
                hcb.Label.String = ['Time (',units,')'];
                hold off; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25 0.25 0.5 0.5]);
                % Plot the other stuff with a nice rainbow and nice axes
                figure; hold on;
                cmap = colormap(jet(size(deltaymatrixselect,2)));
                for i=1:size(deltaymatrixselect,2)
                    plot(xvector(minI:maxI),deltaymatrixselect(2:end,i),'Color',cmap(i,:));
                end
                hline = refline(0,0); hline.Color = [0.5 0.5 0.5];
                set(gca,'FontSize',14)
                xlabel('Wavelength (nm)','FontSize',14,'FontWeight','bold')
                ylabel('\DeltaAbs (OD)','FontSize',14,'FontWeight','bold')
                xlim([xvector(minI),xvector(maxI)]); caxis([0 deltaymatrixselect(1,end)]); hcb=colorbar;
                n = round(deltaymatrixselect(1,end)./6);
                hcb.YTick = 0:n:deltaymatrixselect(1,end);
                hcb.Label.String = ['Time (',units,')'];
                hold off; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25 0.25 0.5 0.5]);
            end
        dooffset = questdlg('Calculate and subtract offset?', ...
        'Offset', ...
        'Yes, please','Yes and plot','No',Struct);
        calc=0;makeplot=0;
        switch dooffset
            case 'Yes, please'
                calc=1; makeplot=0;
            case 'Yes and plot'
                calc=1; makeplot=1;
        end
        if calc == 1
            % Ask for a reduced range and save it
            offsetWL = inputdlg('Define the offset region:','Define offset region',1,{'650 700'});
            offsetWL = str2num(offsetWL{:}); minoffsetWL=offsetWL(1); maxoffsetWL=offsetWL(2);
            name2 = ['DELTA_ABS_' num2str(minoffsetWL) ' to ' num2str(maxoffsetWL) 'nm_OFFSET.dat'];
            deltaymatrixseloffset=[];
            [minJ, ~] = findClosestId2Val(xvector(minI:maxI),minoffsetWL);
            [maxJ, ~] = findClosestId2Val(xvector(minI:maxI),maxoffsetWL);
            offset = mean(deltaymatrixselect(minJ:maxJ,:),1);
            deltaymatrixseloffset = deltaymatrixselect - offset;
            deltaymatrixseloffset(1,:) = deltaymatrixselect(1,:);
            csvwrite(name2,[vertcat(0,xvector(minI:maxI)),deltaymatrixseloffset])
            if makeplot==1
            % Plot the other stuff with a nice rainbow and nice axes
                figure; hold on;
                cmap = colormap(jet(size(deltaymatrixseloffset,2)));
                for i=1:size(deltaymatrixseloffset,2)
                    plot(xvector(minI:maxI),deltaymatrixseloffset(2:end,i),'Color',cmap(i,:));
                end
                hline = refline(0,0); hline.Color = [0.5 0.5 0.5];
                set(gca,'FontSize',14)
                xlabel('Wavelength (nm})','FontSize',14,'FontWeight','bold')
                ylabel('\DeltaAbs (OD)','FontSize',14,'FontWeight','bold')
                xlim([xvector(minI),xvector(maxI)]); caxis([0 deltaymatrixseloffset(1,end)]); hcb=colorbar;
                n = round(deltaymatrixseloffset(1,end)./6);
                hcb.YTick = 0:n:deltaymatrixseloffset(1,end);
                hcb.Label.String = ['Time (',units,')'];
                hold off; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25 0.25 0.5 0.5]);
            end
            tdoffset = questdlg('Show the time evolution of the offset?', ...
            'Offset ticdme evolution', ...
            'Yes','No','No');
            calc=0;makeplot=0;
            switch tdoffset
                case 'Yes'
                figure;
                plot(header,offset,'Color','Red')
                set(gca,'FontSize',14)
                xlabel(['Time (',units,' )'],'FontSize',14,'FontWeight','bold')
                ylabel('\DeltaAbs (OD)','FontSize',14,'FontWeight','bold')
                axis tight
            end
        end
        end
    end
end

function [Ids, DiffVals] = findClosestId2Val(X,Vals)
% Finds the index of the most similar value in an array.
%
% @Input:
%   X (required)
%     A nx1 column vector with the data.
%
%   Vals (required)
%     A 1xk finite row vector with values for which the indices should be found.
%     If 'Vals' contains non-finite values a 'DJM:NonFiniteValues'-warning will
%     be thrown. The according ids will be 1 for -INF, length(X) for INF and NaN
%     for NaN. So check the ids in case of warnings!
%
% @Output:
%   Ids
%     A 1xk vector with the index of the closest values of 'Vals' in 'X'. Ids =
%     1 or Ids = length(X) indicates, that there is no good match to 'Vals' in
%     'X' of some values where infinite.
%
%   DiffVals
%     A 1xk vector with the absolute differences between 'Vals' and the
%     closest value in 'X' (e.g. abs(X(Ids)-Vals)).
%
% @Remarks:
%   - No checking is done on the dimensions of input arguments, so choose your
%     inputs wisely!
%
% @Dependancies
%   none
%
% @Changelog
%   2011-11-13 (DJM): Release.
%   2015-09-29 (DJM): Added support for multiple values.
%   2015-10-28 (DJM): Added some improvements for speed.
%   2016-04-19 (DJM): Added the handling of non-finite values in 'Vals'.
%   2016-08-08 (DJM):
%     - Ids will now be 1 for -INF and length(X) for INF values in 'Vals' as
%       this is closer to the expected behavior. NaNs will still produce NaNs in
%       the ids.
%     - The according warning is been split into a NaN and an INF version. The
%       identifier is still the same for backwards compatibility.

%% Check input
%Assure column vector.
if isrow(X)
	X = X';
end

%Assure row vector.
if iscolumn(Vals)
	Vals = Vals';
end

%Check all-finite.
Ids = NaN(size(Vals));
IsFinite = isfinite(Vals);
if ~all(IsFinite)
	IsInf = isinf(Vals);
	if any(IsInf)
		warning('DJM:NonFiniteValues', ['''Vals'' contained %d infinite ' ...
			'value(s)! The according id(s) will be 1 for -INF or length(X) ' ...
			'for INF!'], nnz(IsInf));
		Ids(Vals == -inf) = 1;
		Ids(Vals ==  inf) = length(X);
	else %NaNs
		warning('DJM:NonFiniteValues', ['''Vals'' contained %d NaN value(s)! ' ...
			'The according id(s) will also be set to NaN!'], nnz(~IsFinite));
	end
end

%% Compute
if any(IsFinite)
	%Get the diffs.
	if isscalar(Vals)
		AbsDiff = X - Vals(IsFinite); %Slightly faster than BSXFUN with scalar.
	else
		AbsDiff = bsxfun(@minus,X,Vals(IsFinite));
	end
	AbsDiff = abs(AbsDiff);

	%Compute the ids.
	[MinDiffVals, Ids(IsFinite)] = min(AbsDiff);
end

%Get the diff values if requested.
if nargout() == 2
	DiffVals = NaN(size(Vals));
	DiffVals(IsFinite) = MinDiffVals;
	if ~all(IsFinite)
		DiffVals(IsInf) = inf; %Since abs we only have positive inf vals.
	end
end