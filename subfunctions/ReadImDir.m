function [filenames] = ReadImDir(directory,ext)
%READIMDIR 
% Reads images in a directory, and returns them as a structure of filenames
% as strings
% [filenames] = ReadImDir(directory,extension)
%
%  Author:  Daniel Buscombe
%           School of Marine Science & Engineering
%           University of Plymouth, Drake Circus, Plymouth, Devon PL48AA, UK
%           daniel.buscombe@plymouth.ac.uk
%  Version: Beta        Revision: 16 April, 2010

direc=dir([directory,filesep,'*.',ext]); %list directory and separate .*ext files
filenames={};   %create a structure of these files
[filenames{1:length(direc),1}] = deal(direc.name); %Deal inputs to outputs!
filenames=sortrows(char(filenames{:})); %Create character array, and sort rows in ascending order
