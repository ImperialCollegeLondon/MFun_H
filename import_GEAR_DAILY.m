

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


function [XX,YY,RAIN] = import_GEAR_DAILY(X_coor,Y_coor,date0,pl)
% Processing GEAR.nc data, for derivation of spatial rainfall pattern.
%
% Input:
%     X_coor
%     Y_coor
%     date0
%     varargin
%
% Output:
%     XX
%     YY
%     RAIN_12
%
%
% Output: monthly spatial pattern of rainfall pattern. 
% RAIN: UNIT: rainfall amount / 365 day
%        
%% save it in different year and different month;

% CHECK no estimate: -2

f_dialog = waitbar(0,'Loading Data ...');

dat = datevec(date0);
year = unique(dat(:,1));
if numel(year)>1
    error('check FUNCTION IMPORT_GEAR_DAILY: input year')
end

dayi = -datenum(datetime(year,1,1,0,0,0))+date0+1;

try
    
    source = ['D:\CEH_GEAR_Daily\CEH_GEAR_daily_GB_',num2str(year),'.nc'];
    finfo = ncinfo(source);
    
    varname = 'rainfall_amount';
    rainfall_amount = ncread(source,varname);
    
    varname = 'x';
    x = ncread(source,varname); % easting-OSGB36 Grid reference
    varname = 'y';
    y = ncread(source,varname); % northing-OSGB36 Grid reference

catch
    fprintf('Check: FUNCTION IMPORT_GEAR_DAILY: DAILY ESTIMATE IS EMPTY');
end

clear year source varname varname mon

%% Define the size;
pause(0.1);

minX = min(X_coor(:));
maxX = max(X_coor(:));
minY = min(Y_coor(:));
maxY = max(Y_coor(:));

indx = find(x<maxX & x>minX);
indy = find(y<maxY & y>minY);

[xx0,yy0] = meshgrid(x,y);

try
    dist = @(x,x0)abs(x-x0);
    indx_min = find(dist(x,minX) == min(dist(x,minX)));
    indx_max = find(dist(x,maxX) == min(dist(x,maxX)));
    indy_min = find(dist(y,minY) == min(dist(y,minY)));
    indy_max = find(dist(y,maxY) == min(dist(y,maxY)));
    
    yi = [indy_max(1),indy_min(1)];
    xi = [indx_min(1),indx_max(1)];
catch
    fprintf('Check: FUNCTION IMPORT_GEAR_DAILY.m');
end

XX = S(yi,xi,xx0);
YY = S(yi,xi,yy0);
RAIN = NaN(size(XX,1),size(XX,2),numel(dayi));
try

for i_d = 1:numel(dayi)
    
    rainfall = rainfall_amount(:,:,dayi(i_d));
    rainfall(rainfall<0) = NaN;
    RAIN(:,:,i_d) = S(yi,xi,rainfall');
    
    waitbar(i_d/numel(dayi),f_dialog,'Processing your data');
    pause(0.01)
    
end

waitbar(1,f_dialog,'Finishing');
pause(1)

catch
    fprintf('Check: FUNCTION IMPORT_GEAR_DAILY.m');
end

if pl
    figure;
    for i = 1:numel(dayi)
        contourf(xx0,yy0,dailyR{i}')
        rectangle('Position',[minX,minY,maxX-minX,maxY-minY]);hold on;
        plot(minX,minY,'wp');
        pause(0.2);
    end
end
delete(f_dialog)

end

function newM = S(Yrange,Xrange,Matrix)

newM = Matrix(Yrange(1):Yrange(2),:);
newM = newM(:,Xrange(1):Xrange(2));

end

%{
%set(gca,'visible','off')
hold on;
for i = 1:length(A)
    a = A(i);
    plot(a.X/unit,a.Y/unit,'--','color','red','linewidth',1);hold on;
end
title('Annual rainfall pettern 2001-2015(GEAR)')
%}

