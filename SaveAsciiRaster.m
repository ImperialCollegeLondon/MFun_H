function SaveAsciiRaster(varname, header,sfPath,sfName);
%   -----  SYNTAX  -----
%   SaveAsciiRaster(varname, header);
%   SaveAsciiRaster(varname);           % in case varname is an xyz matrix
%
%
%   -----  DESCRIPTION  -----
%   This function saves a spatial matrix into an Arc-Info ascii raster. Two file extension '.asc'
%   or '.txt' are supported.
%   FIRST CASE
%       USE: SaveAsciiRaster(varname, header);
%       It requires two inputs: (1) the z-values to be exported ('varname'
%       variable), and (2) the 'header' vector with the spatial information of
%       the grid. 'varname' can be a 1-D vector or a 2-D spatial grid.
%   SECOND CASE
%       USE: SaveAsciiRaster(varname);
%       If an xyz matrix (with [x_coord,y_coord,z_values]) is given as
%       'varname', no 'header' has to be defined, since the function will
%       extract all the required header information from the xyz table. The
%       first row contains the x_coord, y_coord and z_value of the most
%       north-western cell; the last row refers to the most south-eastern
%       pixel. Elements in xyz are sorted column-by-column from the
%       geographical grid
%       (geographical_grid=[1st_col,2nd_col,3rd_col,...,last_col];
%       xyz=[1st_col;2nd_col;3rd_col;...;last_col]).
%       The xy coordinates have to refer to the center of the cells.
%
%
%   -----  INPUT  -----
%  +'varname'(mandatory) : the MatLab matrix to be exported in ascii raster.
%                          It can be passed:
%                          (1) a 2-D z-values matrix of size equal to the spatial grid extent;
%                          (2) a 1-D vector with z-values sorted as [1st_col;2nd_col;...];
%                          (3) an xyz table with characteristics described in "SECOND CASE".
%
%
%  +'header'(facultative): the geospatial reference matrix.
%                          If varname is an xyz matrix then header have not to be given.
%                          The header matrix as the .hdr file in Arc-Info Binary Raster.
%                          The 'header' MatLab variable is created when importing an ascii 
%                          raster with function ImportAsciiRaster (author: Giuliano Langella).
%
%
% Example:
% load(['C:\Users\Yuting Chen\Dropbox (Personal)\MATLAB_CODE_RRMODEL_0319',...
%     '\MATLAB_CODE_RRMODEL\DATA_INPUT\Stratford_TOPO_DATA_200m.mat'],'SN','x_yr','y_yr');
% 
% 
% header = [230; 154; min(x_yr)-100; min(y_yr)-100; 200; -9999];
% sfPath = ['D:\DATA_CAT++\'];
% sfName = ['SN_Stratford.asc'];
% 
% varname = SN;
% SaveAsciiRaster(varname, header,sfPath,sfName);
%
%
%
%
% % Save path
% [FileName, PathName] = uiputfile({'*.asc','Arc-Info ASCII Raster (*.asc)'; ...
%    '*.txt','Text ASCII Raster (*.txt)'; ...
%    '*.*',  'All Files (*.*)'}, 'Save Ascii Grid', pwd);
warning off;
mkdir([sfPath]);
warning on;
PathName = sfPath;
FileName = sfName;

% OPEN file
fid = fopen(strcat(PathName, FileName),'w');
% if I loaded an xyz matrix as 'varname'
if nargin == 1
    cellsize = abs(varname(1,2)-varname(2,2));
    hor = max(varname(:,1)) - min(varname(:,1)) + cellsize;
    ver = max(varname(:,2)) - min(varname(:,2)) + cellsize;
    ncols = ceil(hor/cellsize);
    nrows = ceil(ver/cellsize);
    raster = zeros(nrows,ncols);
    %create the header variable [default no data value is -9999]
    header = [ncols; nrows; min(varname(:,1))-0.5*cellsize; min(varname(:,2))-0.5*cellsize; cellsize; -9999];
    Zvar = reshape(varname(:,3),header(2),header(1));
    varname=[];
    varname=Zvar;
end
% WRITE HEADER
fprintf(fid,'%s','nCols        ');  %1
fprintf(fid,'%12.0f\n', header(1,1));
fprintf(fid,'%s','nRows        ');  %2
fprintf(fid,'%12.0f\n', header(2,1));
fprintf(fid,'%s','xllCorner    ');  %3
fprintf(fid,'%f\n', header(3,1));
fprintf(fid,'%s','yllCorner    ');  %4
fprintf(fid,'%f\n', header(4,1));
fprintf(fid,'%s','CellSize     ');  %5
fprintf(fid,'%f\n', header(5,1));
fprintf(fid,'%s','nodata_value ');  %6
fprintf(fid,'%f\n', header(6,1));
% WRITE MATRIX
%substitute to NaN the NODATA_value written in header
varname(isnan(varname)) = header(6,1);
%start loop
ncols = header(1,1);
nrows = header(2,1);
for CurrRow = 1:nrows;
    % if varname is a vector instead of a 2-D array
    if size(varname,2) == 1;
        fprintf(fid,'% f ',varname( ((CurrRow-1)*ncols + 1) : (CurrRow*ncols) )' );
        fprintf(fid,'%s\n', ' ');
    % if varname is a 2-D array
    else
        fprintf(fid,'%f ',varname(CurrRow,:));
        fprintf(fid,'%s\n', ' ');
    end
end
fclose(fid);
fclose('all');
