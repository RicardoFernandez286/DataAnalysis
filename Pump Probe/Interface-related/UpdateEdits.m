function handles = UpdateEdits(handles)
% Set the edit field labels
handles.editTmin.String         = num2str(handles.plotranges(1));
handles.editTmax.String         = num2str(handles.plotranges(2));
handles.editWLmin.String        = num2str(handles.plotranges(3));
handles.editWLmax.String        = num2str(handles.plotranges(4));
handles.editMaxZ.String         = num2str(handles.plotranges(5));
handles.editNcontours.String    = num2str(handles.plotranges(6));