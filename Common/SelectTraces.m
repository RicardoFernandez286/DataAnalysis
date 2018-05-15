function handles = SelectTraces(handles,doSort)
% Output is a sorted vector (handles.SelTraces) of the selected (X,Y) points via interactive input through the mouse.
%
% Left-click to add points (right-click to remove last)
% Press Return key when done.
%
% Ricardo Fernández-Terán
% v1.1c - 31.03.2018

% Clear the output variable (in case it already existed)
handles.SelTraces = [];

if nargin < 2
    doSort = 1;
end

i=0;
hold on
while 1 % Repeat the loop until Return is pressed
    [xp,yp,button] = ginput(1);
    if isequal(button,1)
        cmap = colormap;
        i=i+1;
        plot(xp,yp,'o','Linewidth',1.5,'Color',cmap(i+1,:))
        handles.SelTraces(i,:)=[xp yp];
    elseif isequal(button,3) && i>0
        % Remove last point when right clicking
        plot(handles.SelTraces(i,1),handles.SelTraces(i,2),...
            'o','Linewidth',1.5,'Color',[1 1 0])
        i=i-1;
    elseif isempty(xp)
            break
    end
end
handles.SelTraces = handles.SelTraces(1:i,:);

if doSort == 1
    handles.SelTraces = sort(handles.SelTraces,1);
end

hold off