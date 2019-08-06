function [dataStruct,aborted] = SelectTraces(dataStruct,doSort,varargin)
% Output is a sorted (or not) vector, dataStruct.SelTraces, of the selected (X,Y) points via interactive input through the mouse.
%
% USAGE:
% Left-click to add points (right-click to remove last)
%
% Press Return key when done.
%
% Press middle button of the mouse to exit
%    (will give an error in the caller routine if nothing was selected)
%
% Ricardo Fernández-Terán
% v2.1a - 09.04.2019

% Clear the output variable (in case it already existed)
dataStruct.SelTraces = [];

warning('off','MATLAB:ginput:FigureDeletionPause');
warning('off','MATLAB:ginput:Interrupted');

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
        dataStruct.SelTraces(i,:)=[xp yp];
    elseif isequal(button,3) && i>0
        % Remove last point when right clicking
        plot(dataStruct.SelTraces(i,1),dataStruct.SelTraces(i,2),'x','Linewidth',1.5,'Color',[1 1 0])
        i=i-1;
    elseif isequal(button,2)
        aborted = 1;
        break
    elseif isempty(xp)
        aborted = 0;
        break
    end
end

dataStruct.SelTraces = dataStruct.SelTraces(1:i,:);

if doSort == 1
    dataStruct.SelTraces = sort(dataStruct.SelTraces,1);
end
hold off
