%function handles = LoadDataIR4(handles,nameformat)

% This function will sort all the files contained in a Lab 4 data directory according to the
% date/time format currently in use for saving the data, then load the data and prepare it for
% further processing.
% 
% INPUTS:
%     handles = structure from DataAnalysis_GUI
%     nameformat = a string containing the current name format (including date and other characters)
%     
% OUTPUTS: handles.$varname. Data are selected according to "sel_datatype"
%     filelist = Cell array containing the filenames of all the valid files to be read.
%                Can be indexed by calling filelist(Scanindices)
%
% Ricardo Fernández-Terán
% v1.0 - 09.01.2017

%% DEBUG
rootdir='D:\Ricardo Data\switchdrive\Ph.D. UZH\RESULTS\Transient absorption\Lab 4 - Pump-probe ATR\180110';
%nameformat=;

% %% READ from handles
% rootdir = handles.rootdir;

%% WRITE to handles
    handles.delays = delays;
    handles.cmprobe = cmprobe;
    handles.rawsignal = rawsignal;
    handles.noise = noise;
    handles.plotranges = plotranges;