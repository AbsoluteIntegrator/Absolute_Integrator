% The intention of this code is to offer as wider range of file types to
% be read by matlab and repackaged ready for passing to another piece of
% code.

% * Note, when using the .mat file format the .mat file must only contain a
% single variable named 'input_file'.


% *************************************************************************


% First the code begins by clearing both the command window and workspace.
clc

% Next a file is selected using the GUI window.
[filename, pathname] = uigetfile( ...
    {'*.*',  'All Files (*.*)'; ...
     '*.dm3','Digital Micrograph Images (*.dm3)'; ...
     '*.mat','MATLAB variables (*.mat)'; ...
     '*.ser','FEI Micrograph Files (*.ser)'; ...
     '*.sti','Simulation Files (*.sti)'; ...
     '*.txt','Text Files (*.txt)'; ...
     '*.bmp;*.png;*.tif;*.jpg;*.gif','Image Files (*.bmp,*.png,*.tif,*.jpg,*.gif)'}, ...
     'File Opener v1.5.2 - Select a file to open...');

% Then its relative filename is converted to an absolute one:
 completename = fullfile(pathname, filename);
 clear pathname
 
% Filenames with more than one '.' charachter are not supported and an
% apropriate error message is generated.
if size(strfind(filename, '.'),1) > 1
    disp( 'Filename with more than one full-stop "." are not valid. Please rename the desired variable.')
end

% The filename of the selected file is scanned for various supported file types.

if size(strfind(filename, '.'),1) == 0 % For files with no file extension generated from simulation.
    file_type = '.bin' ;
elseif strfind(filename, '.bin') > 0
    file_type = '.bin' ;
elseif strfind(filename, '.dm3') > 0
    file_type = '.dm3' ;
elseif strfind(filename, '.gif') > 0
    file_type = '.gif' ;
elseif strfind(filename, '.jpg') > 0
    file_type = '.jpg' ;
elseif strfind(filename, '.JPG') > 0
    file_type = '.jpg' ;
elseif strfind(filename, '.mat') > 0
    file_type = '.mat' ;
elseif strfind(filename, '.txt') > 0
    file_type = '.txt' ;
elseif strfind(filename, '.bmp') > 0
    file_type = '.bmp' ;
elseif strfind(filename, '.tif') > 0
    file_type = '.tif' ;
elseif strfind(filename, '.png') > 0
    file_type = '.png' ;
elseif strfind(filename, '.emi') > 0
    file_type = '.emi' ;
elseif strfind(filename, '.ser') > 0
    file_type = '.ser' ;
elseif strfind(filename, '.sti') > 0
    file_type = '.sti' ;
elseif strfind(filename, '.img') > 0
    file_type = '.img' ;
else
    file_type = 'unknown' ;
end

disp ( horzcat ( 'Selected File is "', filename,'". Detected file type is "',file_type,'".' ) )

if file_type == '.bin'
    
    expression = '_';
    [~,noMatch] = regexp(filename,expression,'match','split') ; % This line splits the filemane into parts delimited by an underscore.
    
    dim = strsplit(noMatch{1,1} , 'x' ) ; % This line splits the first result of the above into two parts delimited by an 'x'.
    
    input_file = fread(fopen(completename), [str2num(dim{1}) str2num(dim{2})] , 'double' , 'ieee-be') ; % This line reads the big-endian data using the dimensions given by 'dim'.
   
elseif file_type == '.dm3'
    input_file = double(ReadDM3(completename))' ;
    
elseif file_type == '.gif'
    input_file = imread (completename) ;
    
elseif file_type == '.jpg'
    input_file = sum(imread (completename),3) ;
    
elseif file_type == '.mat'
    load (completename)
    
elseif file_type == '.txt'
    input_file = importdata (completename) ;
    
elseif file_type == '.bmp'
    input_file = sum(imread (completename),3) ;
    
elseif file_type == '.tif'
    input_file = imread(completename) ;
    input_file = im2double (input_file) ;
    
elseif file_type == '.png'
    [input_file, ~] = imread (completename) ;
    input_file = im2double (input_file) ;
    input_file = sum(input_file,3) ;
    
elseif file_type == '.emi'
    h = msgbox('Sorry ".emi" file-type not presently supported.','File Opener v1.5.5','error') ; 
    
elseif file_type == '.ser'
    serReadIn  = serReader(completename) ;
    input_file = rot90(serReadIn.image)  ;
    clear serReadIn
    
elseif file_type == '.sti'
    DELIMITER = ' ';
    HEADERLINES = 2;
    input_struct = importdata(completename, DELIMITER, HEADERLINES);
    input_file = input_struct.data ;
    clear DELIMITER
    clear HEADERLINES
    clear input_struct
        
elseif file_type == '.img'    
    [input_file,~,~,~] = binread2D(completename) ;
    
end

clear file_type
clear completename