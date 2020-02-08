function IR_SPEC()
% CURRENT VERSION: v2.0
% 04.10.2019 / Ricardo Fernandez-Teran
% Changelog:
% - Updated readParam reference
% - Now it can read Chronoamperometry and CV-type IR-SEC experiments

%% Get the file list
% Get the directory
% read defaultdir, otherwise just show the default dialog
fileexist = exist([fileparts(fileparts(mfilename('fullpath'))) filesep 'GUIoptions.txt'],'file');
if fileexist == 2
    defaultdir = readParam('defaultSPECdir',[fileparts(fileparts(mfilename('fullpath'))) filesep 'GUIoptions.txt']);
    SPECDir = uigetdir(defaultdir);
else
    SPECDir = uigetdir();
end

% If the user selected a directory
if SPECDir ~= 0
ExpType = questdlg('What kind of experiment are we dealing with here?', ...
	'Select IR-SEC type', ...
	'Chronoamperometry','CV','Cancel','CV');
    
% Select only files that are NOT *.dpt, *.txt, *.dat or *.csv
    list = dir(SPECDir); list(1:2,:) = [];
    allfiles = transpose({list.name});
    I = contains(allfiles,["dpt","csv","txt","dat"]);
    files = sort(allfiles(I==0));
    cd(SPECDir);

%% Do something according to the experiment type    
switch ExpType
    case 'Chronoamperometry'       
        % Read the potentials from "Potentials.txt"
        potentials = csvread('Potentials.txt');
        if strcmp(questdlg('Convert to a different potential reference scale?'),'Yes')
            opts.Interpreter = 'tex';
            ref_pot = inputdlg({'Which scale?: ';'E_{1/2} in current scale (mV): '},'Reference potentials',[1,50],{'Fc^{+}/Fc','500'},opts);
            ref_name= ['mV vs ' ref_pot{1}];
            ref_val = str2num(ref_pot{2});
        else
            ref_name= 'mV';
            ref_val = 0;
        end
        potentials   = potentials - ref_val;
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
        [x,id] = sort(x,'ascend');
        ymatrix = vertcat(header,ymatrix(id,:));
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
        catch
            [y,x,~] = ImportOpus(ZeroFile,'RatioAbsorptionChanged');
        end
        x=transpose(x); x=double(x); y=double(y)*1000; %Baseline also in mOD
        [x,id] = sort(x,'ascend');
        baseline = vertcat(0,y(id));
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
            legendtext = [legendtext; {[num2str(potentials(i),'%4.2f') ' mV']}];
        end
        hline = refline(0,0); hline.Color = [0.5 0.5 0.5];
        set(gca,'FontSize',14)
        xlabel('Wavenumbers (cm^{-1})','FontSize',14,'FontWeight','bold')
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
        xlabel('Wavenumbers (cm^{-1})','FontSize',14,'FontWeight','bold')
        ylabel('\DeltaAbs (mOD)','FontSize',14,'FontWeight','bold')
        axis tight; xlim([xvector(minI),xvector(maxI)]); legend(legendtext); legend('boxoff');
        hold off; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25 0.25 0.5 0.5]);
        % Ask the user what to do next
            Struct.Interpreter = 'tex';
            Struct.WindowStyle = 'modal';
            Struct.Default = 'No';
        doOffset = questdlg('Calculate and subtract offset?', ...
            'Offset', ...
            'Yes, please','Yes and plot','No',Struct);
            calc=0;makeplot=0;
            switch doOffset
                case 'Yes, please'
                    calc=1; makeplot=0;
                case 'Yes and plot'
                    calc=1; makeplot=1;
            end
            if calc == 1
                % Ask for a reduced range and save it
                offsetWL = inputdlg('Define the offset region:','Define offset region',1,{'2180 2200'});
                offsetWL = str2num(offsetWL{:}); minoffsetWL=offsetWL(1); maxoffsetWL=offsetWL(2);
                name_off = ['DELTA_ABS_' num2str(minoffsetWL) ' to ' num2str(maxoffsetWL) ' cm-1_OFFSET.dat'];
                deltaymatrixseloffset=[];
                [minJ, ~] = findClosestId2Val(xvector(minI:maxI),minoffsetWL);
                [maxJ, ~] = findClosestId2Val(xvector(minI:maxI),maxoffsetWL);
                offset_avg = mean(deltaymatrixselect(minJ:maxJ,:),1);
                deltaymatrixseloffset = deltaymatrixselect - offset_avg;
                deltaymatrixseloffset(1,:) = deltaymatrixselect(1,:);
                csvwrite(name_off,[vertcat(0,xvector(minI:maxI)),deltaymatrixseloffset])
                if makeplot==1
                % Plot the other stuff with a nice rainbow and nice axes
                    figure; hold on;
                    cmap = colormap(jet(size(deltaymatrixseloffset,2)));
                    for i=1:size(deltaymatrixseloffset,2)
                        plot(xvector(minI:maxI),deltaymatrixseloffset(2:end,i),'Color',cmap(i,:),'LineWidth',2);
                    end
                    hline = refline(0,0); hline.Color = [0.5 0.5 0.5];
                    set(gca,'FontSize',14)
                    xlabel('Wavenumbers (cm^{-1})','FontSize',14,'FontWeight','bold')
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
                    plot(header,offset_avg,'-o','Color','Red')
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
                    ylabel('Wavenumbers (cm^{-1})','FontSize',14,'FontWeight','bold')
                    zlabel('\DeltaAbs (mOD)','FontSize',14,'FontWeight','bold')
                    axis tight
                end
            end
    case 'CV'
        %% Read all files and write everything in a matrix
        xmatrix = [];
        ymatrix = [];
        header  = [];
        time_str= cell(length(files),1);
        for i=1:length(files)
            currfile    = char(files(i));
            IRfile      = [SPECDir filesep currfile];
            try
                [y,x,params] = ImportOpus(IRfile,'RatioAbsorption');
            catch
                [y,x,params] = ImportOpus(IRfile,'RatioAbsorptionChanged');
            end
            xmatrix(:,i)    = x';
            ymatrix(:,i)    = double(y);
            time_str{i}     = strsplit(params.RatioDataAbsorptionChanged.TIM,' ');
            time_num(i,:)   = datenum(time_str{i}(1),'HH:MM:SS.FFF');
        end
        %% Convert relative time into potentials
        %%% 1) Read the saved scan parameters
            delta_t = (time_num - time_num(1))*24*60*60;                                    % s
            Vertex1 = str2double(readParam('Vertex1',[SPECDir filesep 'ScanParams.txt']));  % mV
            Vertex2 = str2double(readParam('Vertex2',[SPECDir filesep 'ScanParams.txt']));  % mV
            Estep   = str2double(readParam('Estep',[SPECDir filesep 'ScanParams.txt']));    % mV
            Estart  = str2double(readParam('Estart',[SPECDir filesep 'ScanParams.txt']));   % mV
            Nscans  = str2double(readParam('Nscans',[SPECDir filesep 'ScanParams.txt']));   % #
            Scanrate= str2double(readParam('Scanrate',[SPECDir filesep 'ScanParams.txt'])); % mV/s
            range   = abs(Vertex1-Vertex2);                                                 % mV
        
        %%% 2) Calculate the time and potential axes
            t1      = (Vertex1-Estart)/Scanrate;
            t2      = (Vertex2-Vertex1)/Scanrate;
            t3      = (Estart-Vertex2)/Scanrate;
            time_ax = linspace(0,abs(t1)+abs(t2)+abs(t3),length(delta_t)*2);
            time_ax = time_ax(1:end-1);
            Volt_t1 = time_ax(time_ax <= abs(t1))*sign(t1)*Scanrate + Estart;
            Volt_t2 = (time_ax(time_ax > abs(t1) & time_ax <= abs(t1)+abs(t2))-abs(t1))*sign(t2)*Scanrate + Vertex1;
            Volt_t3 = (time_ax(time_ax > abs(t1)+abs(t2) & time_ax <= abs(t1)+abs(t2)+abs(t3))-(abs(t1)+abs(t2)))*sign(t3)*Scanrate + Vertex2;
            volt_ax = repmat([Volt_t1 Volt_t2 Volt_t3]',Nscans,1);  %(!)
            dt      = time_ax(2) - time_ax(1);
            time_axV= dt*(1:length(time_ax)*Nscans)-dt;             %(!)%

        %%% 3) Match the delta_t from the experiment to a potential (closest match)
            delta_t = delta_t(delta_t <= max(time_axV));
            pot_idx = findClosestId2Val(time_axV,delta_t);
            pot_t   = volt_ax(pot_idx);
            Nplots  = length(delta_t);
        
        %%% 4) Convert to arbitrary potential reference scale
            if strcmp(questdlg('Convert to a different potential reference scale?'),'Yes')
                opts.Interpreter = 'tex';
                ref_pot = inputdlg({'Which scale?: ';'E_{1/2} in current scale (mV): '},'Reference potentials',[1,50],{'Fc^{+}/Fc','500'},opts);
                ref_name= ['mV vs ' ref_pot{1}];
                ref_val = str2num(ref_pot{2});
            else
                ref_name= 'mV';
                ref_val = 0;
            end
            pot_t   = pot_t - ref_val;
            
        %%% 5) Discard the columns corresponding to data after the CV finished
            xmatrix = xmatrix(:,1:Nplots);      % in cm-1
            xvector = mean(xmatrix,2);                   % in cm-1
            ymatrix = ymatrix(:,1:Nplots)*1000; % in mOD

            csvwrite('ALLDATA.dat',[[0, pot_t'];[xvector,ymatrix]]);
        
        %%% 6) Plot the raw data (if requested)
            if strcmp(questdlg('Plot raw data?','Raw data','Yes, please','No','No'),'Yes, please')
                fh  = figure;
                fh.Color = 'w';
                ax  = axes('parent',fh);
                cmap= colormap(jet(Nplots));
                hold(ax,'on');
                for i=1:Nplots
                    plot(ax,xvector,ymatrix(:,i),'-','Color',cmap(i,:),'DisplayName',num2str(pot_t(i),'%2.f'));
                end
                hold(ax,'off');
                axis(ax,'tight');
                xlabel(ax,'Wavenumbers (cm^{-1})','FontWeight','bold','FontSize',14);
                ylabel(ax,'Abs (mOD)','FontWeight','bold','FontSize',14);
                box(ax,'on');
                ax.FontSize = 14;
                legend(ax,'show');
                legend(ax,'Location','bestoutside');
                legend(ax,'boxoff');
            end
            
        %%% 6) Inform the user
            uiwait(msgbox('Conversion completed! Now you will be asked to select the baseline file', 'Success!','warn'));
            
        %% Plot the CV difference data (with and without offset correction)
        %%% 1) Ask for the open potential (or 0 V) scan
            [ZeroPotFile,ZeroPath,~] = uigetfile({'*.?*','Bruker OPUS files (*.#)'},'Select the baseline spectrum to load...');
            if ZeroPotFile == 0
                return
            end
            % Use this as the "zero potential" baseline
            try
                [y,x,~] = ImportOpus([ZeroPath ZeroPotFile],'RatioAbsorption');
            catch
                [y,x,~] = ImportOpus([ZeroPath ZeroPotFile],'RatioAbsorptionChanged');
            end
            baseline_x  = x';
            baseline_y  = double(y)*1000; %Baseline also in mOD

            if length(baseline_x) ~= length(xvector) || norm(baseline_x - xvector) >= 1e-5
                errordlg('Inconsistent baseline spectrum!','Error');
                return
            end
            
        %%% 2) Calculate the DeltaY matrix, plot it and save it
            deltaYmatrix = ymatrix - baseline_y;
            csvwrite('DELTA_ABS.dat',[[0, pot_t'];[mean(xmatrix,2),deltaYmatrix]]);
            
        %%% 4) Ask for a reduced range to save (and plot) the data with and without the offset
            % Get
            WL          = inputdlg('Define the spectral window to save:','Define spectral window',1,{'1775 2200'});
            WL          = str2num(WL{:}); minWL=WL(1); maxWL=WL(2);
            name_sel    = ['DELTA_ABS_' num2str(minWL) ' to ' num2str(maxWL) ' cm-1.dat'];
            [minI, ~]   = findClosestId2Val(xvector,minWL);
            [maxI, ~]   = findClosestId2Val(xvector,maxWL);
            deltaYmatrix_sel    = deltaYmatrix(minI:maxI,:);
            xvector_sel         = xvector(minI:maxI);
            csvwrite(name_sel,[[0, pot_t'];[xvector_sel,deltaYmatrix_sel]]);
          
        %%% 5) Subtract the offset (if wanted) and save it
            doOffset            = questdlg('Calculate and subtract offset?','Offset','Yes, please','Yes and plot','No','No');
            calc                = 0;
            makeOffsetPlot      = 0;
            
            switch doOffset
                case 'Yes, please'
                    calc=1; makeOffsetPlot=0;
                case 'Yes and plot'
                    calc=1; makeOffsetPlot=1;
            end
            
            if calc == 1
                % Ask for an offset range, calculate the average and subtract it
                offsetWL    = inputdlg('Define the offset region:','Define offset region',1,{'2180 2200'});
                offsetWL    = str2num(offsetWL{:});
                minoffsetWL = offsetWL(1);
                maxoffsetWL = offsetWL(2);
                name_off    = ['DELTA_ABS_' num2str(minoffsetWL) ' to ' num2str(maxoffsetWL) ' cm-1_OFFSET.dat'];
                [minJ, ~]   = findClosestId2Val(xvector,minoffsetWL);
                [maxJ, ~]   = findClosestId2Val(xvector,maxoffsetWL);
                offset_avg  = mean(deltaYmatrix(minJ:maxJ,:),1);
                deltaYmatrix_selOff = deltaYmatrix_sel - offset_avg;
                csvwrite(name_off,[[0, pot_t'];[xvector_sel,deltaYmatrix_selOff]]);
            else
                offset_avg  = 0;
            end
            
            plotDA = strcmp(questdlg('Plot full DeltaAbs data?','Full DeltaAbs data','Yes, please','No','No'),'Yes, please');
            
            if plotDA == 1
                fh  = figure;
                fh.Color = 'w';
                ax  = axes('parent',fh);
                cmap= colormap(jet(Nplots));
                hold(ax,'on');
                for i=1:Nplots
                    plot(ax,xvector,deltaYmatrix(:,i)-offset_avg,'-','Color',cmap(i,:),'DisplayName',num2str(pot_t(i),'%2.f'));
                end
                hold(ax,'off');
                axis(ax,'tight');
                set(fh, 'Units', 'Normalized', 'OuterPosition', [0.25 0.25 0.5 0.5]);
                xlabel(ax,'Wavenumbers (cm^{-1})','FontWeight','bold','FontSize',14);
                ylabel(ax,'Abs (mOD)','FontWeight','bold','FontSize',14);
                box(ax,'on');
                ax.FontSize = 14;
                legend(ax,'show');
                legend(ax,'Location','bestoutside');
                legend(ax,'boxoff');
            end
            
       %% Plot the stuff 
       %%% 1) Plot the DeltaAbs data
%             if plotDA == 1
%                 fh  = figure;
%                 fh.Color = 'w';
%                 ax  = axes('parent',fh);
%                 cmap= colormap(jet(Nplots));
%                 hold(ax,'on');
%                 for i=1:Nplots
%                     plot(ax,xvector_sel,deltaYmatrix_sel(:,i),'-','Color',cmap(i,:),'DisplayName',num2str(pot_t(i),'%2.f'));
%                 end
%                 hold(ax,'off');
%                 axis(ax,'tight');
%                 set(fh, 'Units', 'Normalized', 'OuterPosition', [0.25 0.25 0.5 0.5]);
%                 xlabel(ax,'Wavenumbers (cm^{-1})','FontWeight','bold','FontSize',14);
%                 ylabel(ax,'\DeltaAbs (mOD)','FontWeight','bold','FontSize',14);
%                 box(ax,'on');
%                 ax.FontSize = 14;
%                 legend(ax,'show');
%                 legend(ax,'Location','bestoutside');
%                 legend(ax,'boxoff');
%             end
            
       %%% 2) Plot the DeltaAbs offset subtracted data
            if makeOffsetPlot == 1
                fh  = figure;
                fh.Color = 'w';
                ax  = axes('parent',fh);
                cmap= colormap(jet(Nplots));
                hold(ax,'on');
                for i=1:Nplots
                    if mod(i,4) == 0
                        plot(ax,xvector_sel,deltaYmatrix_selOff(:,i),'-','Color',cmap(i,:),'DisplayName',num2str(pot_t(i),'%2.f'));
                    else
                        plot(ax,xvector_sel,deltaYmatrix_selOff(:,i),'-','Color',cmap(i,:),'DisplayName',num2str(pot_t(i),'%2.f'),'HandleVisibility','off');
                    end
                end
                hold(ax,'off');
                axis(ax,'tight');
                set(fh, 'Units', 'Normalized', 'OuterPosition', [0.25 0.25 0.5 0.5]);
                xlabel(ax,'Wavenumbers (cm^{-1})','FontWeight','bold','FontSize',14);
                ylabel(ax,'\DeltaAbs (mOD)','FontWeight','bold','FontSize',14);
                box(ax,'on');
                ax.FontSize = 14;
                legend(ax,'show');
                legend(ax,'Location','bestoutside');
                legend(ax,'boxoff');
            end
            
        %% Make a nice contour plot of the offset-subtracted data (if exists), otherwise raw deltaAbs data
        %%% Options
            Ncontours           = 40;
            plot_skiplevels     = 2;
            PlotContourLines    = 1;
            LineColor           = 0.4*[1 1 1];
            ContourLineStyle    = '-';

        %%% Plot
            if calc == 1
                plotwhat    = deltaYmatrix_selOff;
            else
                plotwhat    = deltaYmatrix_sel;
            end

            fh          = figure;
            fh.Color    = 'w';
            ax          = axes('parent',fh);
            min_cut     = -max(abs(plotwhat(:)));
            max_cut     = max(abs(plotwhat(:)));
            cmap        = b2r(min_cut,max_cut,Ncontours,2);

            if plot_skiplevels > 1
                negCont         = decim(linspace(min_cut,0,round((Ncontours+1)/2)),plot_skiplevels,'max');
                posCont         = decim(linspace(0,max_cut,round((Ncontours+1)/2)),plot_skiplevels,'max');
                contourlines    = [negCont; posCont];
                contoursolid    = linspace(min_cut,max_cut,Ncontours+1);
                contoursolid    = contoursolid(1:end-1);
            else
                negCont         = linspace(min_cut,0,round((Ncontours+1)/2))';
                posCont         = linspace(0,max_cut,round((Ncontours+1)/2))';
                contourlines    = [negCont; posCont];
                contoursolid    = linspace(min_cut,max_cut,Ncontours+1);
                contoursolid    = contoursolid(1:end-1);
            end

            contourf(ax,1:1:Nplots,xvector_sel,plotwhat,contoursolid,'LineColor','flat','LineStyle',ContourLineStyle,'LineWidth',0.05);
            if PlotContourLines == 1
                hold(ax,'on')
                contour(ax,1:1:Nplots,xvector_sel,plotwhat,contourlines,'LineColor',LineColor,'LineStyle',ContourLineStyle,'LineWidth',0.05);
                hold(ax,'off')
            end

            colormap(ax,cmap);
            shading(ax,'flat');
            caxis(ax,[-max(abs(plotwhat(:))) max(abs(plotwhat(:)))]);
            axis(ax,'tight');
            xlabel(ax,['Potential (' ref_name ')'],'FontWeight','bold','FontSize',14);
            ylabel(ax,'Wavenumbers (cm^{-1})','FontWeight','bold','FontSize',14);
            xticks(ax,linspace(1,length(pot_t),10));
            xticklabels(ax,num2str(pot_t(round(linspace(1,length(pot_t),10))),'%2.f'));
            colorbar(ax);
            box(ax,'on');
            ax.FontSize = 14;
            ax.Layer    = 'top';
            set(fh, 'Units', 'Normalized', 'OuterPosition', [0.25 0.25 0.5 0.5]);
end
else
    return
end