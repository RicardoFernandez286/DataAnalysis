function handles = UpdateEdits(handles)
% Read from handles
minwl       = num2str(handles.plotranges(3));
maxwl       = num2str(handles.plotranges(4));
Zminmax     = num2str(handles.plotranges(5));
Ncontours   = num2str(handles.plotranges(6));

% Set the edit field labels
set(handles.editWLmin,'String',minwl)
set(handles.editWLmax,'String',maxwl)
set(handles.editMaxZ,'String',Zminmax)
set(handles.editNcontours,'String',Ncontours)