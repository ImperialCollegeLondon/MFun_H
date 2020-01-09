

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORRECT RADAR USING CEH-GEAR DAILY DATA
%
% Sourse File: 'I:\CEH_GEAR_Daily'
% Author: Yuting Chen
% Imperial College London
% Explaination: ..............
% yuting.chen17@imperial.ac.uk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [XX,YY,RAIN] = import_GEAR_DAILY(X_coor, Y_coor, date0, pl)
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
%
% Output:
%     XX
%     YY
%     RAIN_12
% Output: monthly spatial pattern of rainfall
% RAIN: UNIT: rainfall amount / 365 day
% 
%
% Example:
% X_coor = 299000:1000:304000;
% Y_coor = 636000:1000:651000;
% [X_coor,Y_coor] = meshgrid(X_coor,Y_coor);
% 
% date0 = datenum(datetime(2010,1,1):datetime(2010,12,31));
% [XX,YY,RAIN] = import_GEAR_DAILY(X_coor,Y_coor,date0,[]);
%
%
% %%
% YEARRANGE = [2011:2015];
% options = weboptions('username','yutingchen0604@hotmail.com','password','AaBb14207','Timeout',Inf);
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

dat = datevec(date0);
year = unique(dat(:,1));
if numel(year)>1
    error('check FUNCTION IMPORT_GEAR_DAILY: input year')
end

dayi = round(-datenum(datetime(year,1,1,0,0,0))+date0+1);

%D:\CEH_GEAR_Daily
source = ['K:\GEAR\CEH_GEAR_daily_GB_',num2str(year),'.nc'];
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
ti = [dayi(1),dayi(end)];

varname = 'rainfall_amount';
if yi(2)<yi(1)
    RAIN=ncread(source,varname,...
        [xi(1),yi(1),ti(1)],...
        [xi(2)-xi(1)+1,yi(1)-yi(2)+1,ti(2)-ti(1)+1]);
    RAIN = RAIN(:,end:-1:1,:);
else
    RAIN=ncread(source,varname,...
        [xi(1),yi(1),ti(1)],...
        [xi(2)-xi(1)+1,yi(2)-yi(1)+1,ti(2)-ti(1)+1]);
end
XX = X_coor;
YY = Y_coor;

end
