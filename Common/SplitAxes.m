function [axlink, ax_left]=SplitAxes(ax_handle,fig_handle)

% Split an axes in two, the first is around zero and the second is a log
% plot.
% 
% 2015-03-31, Jonas Petersson, peterssonjonas@yahoo.se
% Inspired by Jens Föhlinger

%% Get current axes and figure if they are not input

if not(exist('ax_handle','var')) || isempty(ax_handle)
    ax_handle=gca;
end
if not(exist('fig_handle','var')) || isempty(fig_handle)
    fig_handle=gcf;
end


%% Default Settings (set in plotsettings.txt)

% Time to split axis, and start of the left axis
t_split=str2num(readParam('AxisBreak','C:\GUIOptions.txt')); % e.g. [-0.5 0.5]
if isempty(t_split)
    t_split=get(ax_handle,'XLim');
    t_split(2)=t_split(1)+0.1*(t_split(2)-t_split(1));
end

%Width of the first axis relative the total width,
%and distance between axis
ax_ratio=str2num(readParam('BreakRatio','C:\GUIOptions.txt')); % e.g. [0.25 0.01]
if isempty(ax_ratio)
    ax_ratio=[0.25 0.01];
end

%% Start splitting!

%x-lim of original axes
t=get(ax_handle,'XLim');
%x-lim of the new axes
tlim_left=t_split;
tlim_right=[t_split(2) t(2)];

pos=get(ax_handle,'Position'); %position of original axis
ax_left=copyobj(ax_handle,fig_handle);

%The original axes is moved to the right
pos_right=pos;
pos_right(1)=pos(1)+pos(3)*ax_ratio(1);
pos_right(3)=pos(3)-pos(3)*ax_ratio(1);
set(ax_handle,'Position',pos_right,'XLim',tlim_right,'YTickLabel',[])
set(ax_handle,'XScale',readParam('LinLogScale','C:\GUIOptions.txt'))
set(get(ax_handle,'YLabel'),'String',[])
% set(a1,'XTickLabel',num2str(get(a1,'XTick')'))
yl_right=get(ax_handle,'YLim');

%Duplicated axes to the left
pos_left=pos;
pos_left(3)=pos(3)*ax_ratio(1)-ax_ratio(2)*pos(3);
set(ax_left,'Position',pos_left,'XLim',tlim_left,'XScale','linear')
set(get(ax_left,'XLabel'),'String',[])
set(get(ax_left,'Title'),'String',[])
yl_left=get(ax_left,'YLim');

% link the -y-axes
linkaxes([ax_handle ax_left],'y')
axlink = linkprop([ax_handle,ax_left],'YLim'); %(for multiple linking)

set(ax_handle,'YLim',[min([yl_right yl_left]) max([yl_right yl_left])])
% axes(ax_handle); %bring the right axis (original axis) to focus

