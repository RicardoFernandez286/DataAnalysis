function [yData,xData,params] = ImportOpus(path,strDataType)
%Imports the data from Bruker Opus files into Matlab
%INPUT: Takes single input parameter, 'path' - the full path to the file
%i.e. 'C:\Data\OpusFiles\file.0'
%
%ADDITIONAL INPUT:
%strDataType is a string describing the required data
%type to be loaded.  Use of this parameter removes the requirement of user
%input during the loading process making loading multiple files in a loop
%more simple.
%
%Example data type strings:
%'SampleInterferogram'
%'SampleSpectrum'
%'RatioAbsorption'
%'ReferenceSpectrum'
%
%OUTPUT
%yData - the spectral data
%xData - the wavenumber data (or equivalent)
%params -  struct of structs containing all parameters stored in the
%opus files
%
%This function was written by Jacob Filik and is free to use and
%distribute, if you have any comments or suggestions, please contact:
%
%jacob.filik@diamond.ac.uk
%jacob.filik@gmail.com