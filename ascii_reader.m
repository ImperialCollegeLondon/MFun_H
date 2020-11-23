function [OUT, ncols, nrows, xllcorner, yllcorner, cellsize, nodata] = ascii_reader (filename,varargin)

% [OUT, ncols, nrows, xllcorner, yllcorner, cellsize, nodata] = ascii_reader (filename)

% This function reads arc .asc files
% requires filename string as input

% @Yuting

%% read an ascii file
fin = fopen(filename,'r');
if fin == -1
    error('No file');
end
A = fscanf(fin,'%s',1); ncols = fscanf(fin,'%f',1);           %#ok<NASGU>
A = fscanf(fin,'%s',1); nrows = fscanf(fin,'%f',1);           %#ok<NASGU>
A = fscanf(fin,'%s',1); xllcorner = fscanf(fin,'%f',1);       %#ok<NASGU>
A = fscanf(fin,'%s',1); yllcorner = fscanf(fin,'%f',1);       %#ok<NASGU>
A = fscanf(fin,'%s',1); cellsize = fscanf(fin,'%f',1);        %#ok<NASGU>

if isempty(varargin) || varargin{1} == 6
    A = fscanf(fin,'%s',1); nodata = fscanf(fin,'%f',1);          %#ok<NASGU>
else
    nodata = [];
end

% OUT = fscanf(fin,'%f',[ncols, nrows]);
OUT = cell2mat(textscan(fin,repmat('%f',1,ncols))); % quicker

% OUT = OUT';

fclose('all');
