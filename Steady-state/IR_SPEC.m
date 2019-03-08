function IR_SPEC()
% CURRENT VERSION v1.2
% 04.01.2018 / Ricardo Fernández-Terán
% Changelog:
% - Updated readParam reference

% Get the directory
% read defaultdir, otherwise just show the default dialog
fileexist = exist('C:\GUIoptions.txt');
if fileexist == 2
defaultdir = readParam('defaultSPECdir','C:\GUIoptions.txt');
SPECDir = uigetdir(defaultdir);
else
SPECDir = uigetdir();
end
% If the user selected a directory
if SPECDir ~= 0
    % Select only files that are NOT *.dpt, *.txt, *.dat or *.csv
    list = dir(SPECDir); list(1:2,:) = [];
    allfiles = transpose({list.name});
    I = contains(allfiles,["dpt","csv","txt","dat"]);
    files = sort(allfiles(I==0));
    cd(SPECDir);
        % Read the potentials from "Potentials.txt"
        potentials = csvread('Potentials.txt');
        % Read all files and write everything in a matrix
        ymatrix = []; xvector = []; header = [];
        for i=1:length(files)
            currfile = char(files(i));
            IRfile = [SPECDir,filesep,currfile];
            try
                [y,x,~] = ImportOpus(IRfile,'RatioAbsorption');
            catch err
                [y,x,~] = ImportOpus(IRfile,'RatioAbsorptionChanged');
            end
            x=transpose(x); x=double(x); y=double(y);
            newname = [currfile,'.dat'];
            ymatrix = [ymatrix,y];
        end
        header = transpose(potentials);
        ymatrix = vertcat(header,ymatrix);
        xvector = vertcat(0,x);
        csvwrite('ALLDATA.dat',[xvector,ymatrix])
        uiwait(msgbox('Conversion completed! Now you will be asked to select the baseline file', 'Success!','warn'));
        % Convert to mOD
        ymatrix =  ymatrix * 1000;
        % Ask for the open potential (or 0 V) scan
        [ZeroPotFile,ZeroPath,~] = uigetfile({'*.?*','Bruker OPUS files (*.#)'},'Select the baseline spectrum to load...');
        % Use this as the "zero potential" baseline
        ZeroFile = [ZeroPath ZeroPotFile];
        try
            [y,x,~] = ImportOpus(ZeroFile,'RatioAbsorption');
        catch err
            [y,x,~] = ImportOpus(ZeroFile,'RatioAbsorptionChanged');
        end
        x=transpose(x); x=double(x); y=double(y)*1000; %Baseline also in mOD
        baseline = vertcat(0,y);
        deltaymatrix = [];
        for i=1:size(ymatrix,2) % # of columns of ymatrix
            deltaymatrix = [deltaymatrix,ymatrix(:,i)-baseline];
        end
        deltaymatrix = vertcat(header,deltaymatrix(2:end,:));
        % Save it
        csvwrite('DELTA_ABS.dat',[xvector,deltaymatrix])
        % Ask for a reduced range and save it
        WL = inputdlg('Define the spectral window to save:','Define spectral window',1,{'1775 2200'});
        WL = str2num(WL{:}); minWL=WL(1); maxWL=WL(2);
        name = ['DELTA_ABS_' num2str(minWL) ' to ' num2str(maxWL) ' cm-1.dat'];
        deltaymatrixselect=[];
        [minI, ~] = findClosestId2Val(xvector,minWL);
        [maxI, ~] = findClosestId2Val(xvector,maxWL);
        deltaymatrixselect = vertcat(header,deltaymatrix(minI:maxI,:));
        csvwrite(name,[vertcat(0,xvector(minI:maxI)),deltaymatrixselect])

        % Plot the stuff with a rainbow and nice axes
        figure; hold on;
        cmap = colormap(jet(size(deltaymatrix,2)));
        legendtext = {};
        for i=1:size(deltaymatrix,2)
            plot(xvector(2:end),deltaymatrix(2:end,i),'Color',cmap(i,:),'LineWidth',2);
            legendtext = [legendtext; {[num2str(potentials(i)) ' mV']}];
        end
        hline = refline(0,0); hline.Color = [0.5 0.5 0.5];
        set(gca,'FontSize',14)
        xlabel('Wavenumbers (cm^{�1})','FontSize',14,'FontWeight','bold')
        ylabel('\DeltaAbs (mOD)','FontSize',14,'FontWeight','bold')
        xlim([xvector(2),xvector(end)]); legend(legendtext); legend('boxoff');
        hold off; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25 0.25 0.5 0.5]);

        % Plot the other stuff with a nice rainbow and nice axes
        figure; hold on;
        cmap = colormap(jet(size(deltaymatrixselect,2)));
        for i=1:size(deltaymatrixselect,2)
            plot(xvector(minI:maxI),deltaymatrixselect(2:end,i),'Color',cmap(i,:),'LineWidth',2);
        end
        hline = refline(0,0); hline.Color = [0.5 0.5 0.5];
        set(gca,'FontSize',14)
        xlabel('Wavenumbers (cm^{�1})','FontSize',14,'FontWeight','bold')
        ylabel('\DeltaAbs (mOD)','FontSize',14,'FontWeight','bold')
        axis tight; xlim([xvector(minI),xvector(maxI)]); legend(legendtext); legend('boxoff');
        hold off; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25 0.25 0.5 0.5]);
        % Ask the user what to do next
            Struct.Interpreter = 'tex';
            Struct.WindowStyle = 'modal';
            Struct.Default = 'No';
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
                offsetWL = inputdlg('Define the offset region:','Define offset region',1,{'2180 2200'});
                offsetWL = str2num(offsetWL{:}); minoffsetWL=offsetWL(1); maxoffsetWL=offsetWL(2);
                name2 = ['DELTA_ABS_' num2str(minoffsetWL) ' to ' num2str(maxoffsetWL) ' cm-1_OFFSET.dat'];
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
                        plot(xvector(minI:maxI),deltaymatrixseloffset(2:end,i),'Color',cmap(i,:),'LineWidth',2);
                    end
                    hline = refline(0,0); hline.Color = [0.5 0.5 0.5];
                    set(gca,'FontSize',14)
                    xlabel('Wavenumbers (cm^{�1})','FontSize',14,'FontWeight','bold')
                    ylabel('\DeltaAbs (mOD)','FontSize',14,'FontWeight','bold')
                    axis tight; xlim([xvector(minI),xvector(maxI)]); legend(legendtext); legend('boxoff');
                    hold off; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25 0.25 0.5 0.5]);
                end
                tdoffset = questdlg('Show the evolution of the offset?', ...
                'Offset evolution', ...
                'Yes','No','No');
                switch tdoffset
                    case 'Yes'
                    figure;
                    plot(header,offset,'-o','Color','Red')
                    set(gca,'FontSize',14)
                    xlabel(['Potential (mV)'],'FontSize',14,'FontWeight','bold')
                    ylabel('\DeltaAbs (mOD)','FontSize',14,'FontWeight','bold')
                    axis tight
                end
                plot3D = questdlg('Make a 3D plot?', ...
                '3D plot', ...
                'Yes','No','No');
                switch plot3D
                    case 'Yes'
                    figure;
                    s=surf(deltaymatrixseloffset(1,:),xvector(minI:maxI),deltaymatrixseloffset(2:end,:));
                    s.EdgeColor = 'interp'; s.FaceAlpha=0.6; s.MeshStyle='column';
                    set(gca,'FontSize',14)
                    xlabel(['Potential (mV)'],'FontSize',14,'FontWeight','bold')
                    ylabel('Wavenumbers (cm^{�1})','FontSize',14,'FontWeight','bold')
                    zlabel('\DeltaAbs (mOD)','FontSize',14,'FontWeight','bold')
                    axis tight
                end
            end
end
end