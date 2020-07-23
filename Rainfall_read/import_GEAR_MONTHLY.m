

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORRECT RADAR USING CEH-GEAR MONTHLY DATA
%
% Sourse File: 'K:\GEAR-month'
% Author: Yuting Chen
% Imperial College London
% Explaination: ..............
% yuting.chen17@imperial.ac.uk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [XX,YY,RAIN,DIST] = import_GEAR_MONTHLY(X_coor, Y_coor, year, month, pl)
% Processing GEAR.nc data, for derivation of spatial rainfall pattern.
%
% Input:
%     X_coor
%     Y_coor
%     date0
%     varargin
%
% The input have to be the resolution of 1Km
% X_coor, Y_coor: meshgrid format, with the same distance (1km), can be
% increasig/descending
% The input length has to be within one year (only read one file in this function)
% Also Be careful about the different Coor order.
%
% Output:
%     XX:
%     YY:
%     RAIN:
%     DIST:minimum distance to rain gauges.
% Output: monthly spatial pattern of rainfall
% RAIN: UNIT: rainfall amount / 365 day
%       format: 2d # [Northing, Easting]
% 
%
% Example:
% X_coor = 299000:1000:304000;
% Y_coor = 636000:1000:651000;
% [X_coor,Y_coor] = meshgrid(X_coor,Y_coor);
% 
% [XX,YY,RAIN] = import_GEAR_MONTHLY(X_coor,Y_coor,2010,12,[]);
%
%
% %%
% YEARRANGE = [2011:2015];
% options = weboptions('username',getYutingEmail(),'password',getYutingEmail('PAssword'),'Timeout',Inf);
% weblink = 'https://catalogue.ceh.ac.uk/datastore/eidchub/ee9ab43d-a4fe-4e73-afd5-cd4fc4c82556/GB/daily/';
% for year = YEARRANGE
%
%     try
%         OFN = websave(['K:\GEAR\CEH_GEAR_daily_GB_',num2str(year),'.nc'],...
%             [weblink,'CEH_GEAR_daily_GB_',num2str(year),'.nc'],options);
%         disp(['Finished:',num2str(year)]);
%     catch
%         disp('111');
%     end
%
% end
%
%
%% save it in different year and different month;

arguments
    X_coor (:,:) double
    Y_coor (:,:) double
    year (1,1) double
    month (1,1) double
    pl (1,1) logical = false
end


source = ['K:\GEAR-month\CEH_GEAR_monthly_GB_',num2str(year),'.nc'];
finfo = ncinfo(source);

varname = 'x';
x = ncread(source,varname); % easting-OSGB36 Grid reference
varname = 'y';
y = ncread(source,varname); % northing-OSGB36 Grid reference

%% Define the size;

X_coor = round(X_coor/1000)*1000;
Y_coor = round(Y_coor/1000)*1000;

X1 = X_coor(1,1);
X2 = X_coor(end,end);
Y1 = Y_coor(1,1);
Y2 = Y_coor(end,end);

dist = @(x,x0)abs(x-x0);
get_ind = @(X0,x)find(dist(x,X0) == min(dist(x,X0)));
xi = [get_ind(X1,x),get_ind(X2,x)];
yi = [get_ind(Y1,y),get_ind(Y2,y)];


varname = 'rainfall_amount';
if yi(2)<yi(1)
    RAIN=ncread(source,varname,...
        [xi(1),yi(2),month],...
        [xi(2)-xi(1)+1,yi(1)-yi(2)+1,1]);
    RAIN = RAIN(:,end:-1:1,:);
else
    RAIN=ncread(source,varname,...
        [xi(1),yi(1),month],...
        [xi(2)-xi(1)+1,yi(2)-yi(1)+1,1]);
end

varname = 'min_dist';
if yi(2)<yi(1)
    DIST=ncread(source,varname,...
        [xi(1),yi(2),month],...
        [xi(2)-xi(1)+1,yi(1)-yi(2)+1,1]);
    DIST = DIST(:,end:-1:1,:);
else
    DIST=ncread(source,varname,...
        [xi(1),yi(1),month],...
        [xi(2)-xi(1)+1,yi(2)-yi(1)+1,1]);
end

XX = X_coor;
YY = Y_coor;
RAIN = permute(RAIN,[2,1]);
DIST = permute(DIST,[2,1]);

end



