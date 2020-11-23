
% This file saves all info of regions.
% resolution: 2.2km
function REGIONS = REGIONS_info()

REGIONS = struct;
%% STREAP CALIBRATION RELATED
REGIONS.Cleehill = struct('Name','Cleehill','minE',303,'minN',222,'dimE',110,'dimN',110,'dx',1);


%% UKCP EVALUTION RELATED
% areas used to compare his ukcp & fut ukcp
REGIONS.SCO = struct('Name','bigSCO','minE',50,'minN',550,'dimE',200,'dimN',200,'dx',2.2);
REGIONS.WAL = struct('Name','bigWAL','minE',100,'minN',0,'dimE',140,'dimN',200,'dx',2.2);
REGIONS.EUK = struct('Name','bigEUK','minE',450,'minN',0,'dimE',150,'dimN',150,'dx',2.2);

% areas used to compare radar & his ukcp
REGIONS.SWestuk = struct('Name','SWestuk','minE',221,'minN',41,'dimE',50,'dimN',50,'dx',2.2);
REGIONS.Westuk = struct('Name','Westuk','minE',219.8,'minN',225.8,'dimE',50,'dimN',50,'dx',2.2);
REGIONS.Scotland = struct('Name','Scotland','minE',206.7,'minN',627.6,'dimE',50,'dimN',50,'dx',2.2);
REGIONS.London = struct('Name','London','minE',475,'minN',125.3,'dimE',50,'dimN',50,'dx',2.2);

% .UK is defined as overlapped area of UKCP and Radar.
REGIONS.UK = struct('Name','UK','minE',-280,'minN',-226,'dimE',500,'dimN',640,'dx',2.2);

% region info for Big Birm (all related WCs are included.)
REGIONS.Birm = struct('Name','Birm','minE',390,'minN',270,'dimE',30,'dimN',40,'dx',1);
% region info for East Anglia (one area picked up, no specific criteria)
REGIONS.EAng = struct('Name','EastAnglia','minE',530,'minN',220,'dimE',50,'dimN',50,'dx',2.2);
% BIGGER AREA
% here UK region is the overlapped area of UKCP and Radar.
REGIONS.smallUK = struct('Name','UK_','minE',0,'minN',0,'dimE',315,'dimN',565,'dx',2.2);
% only south of UK
REGIONS.SUK = struct('Name','UK_','minE',0,'minN',0,'dimE',315,'dimN',180,'dx',2.2);
end



%
% plot(borderE/1000,borderN/1000,'k-','linewidth',0.5)
% hold on;rectangle('Position',[50,550,440,440])
% hold on;rectangle('Position',[100,0,308,440]);
% hold on;rectangle('Position',[450,0,330,330]);


