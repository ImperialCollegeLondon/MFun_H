function [lengthi,stormsi] = seperateRadarEvents( radarWAR , thresholdWAR , sepHour, thresholdTime, radarReso)
% -------------------------------------------------------------------------
% seperateRadarEvents Define wet event base on the radar data
%
% Inputs:
%   radarWAR - wer area ratio [-]
%   thresholdWAR - wet area threshold to consider as wet event [-] (recommend? 0.02, used in STREAP)
%   sepHour - minimum seprate hour (recommend: 2, used in STREAP)
%   threshold time - minimum time of wet event lentgh [min] (recommend: ??)
%   radarReso - resolution of radar data [min] 
%               (default radar reso is: 5min, please not coarser than 60 minutes)
% Output:
%   lengthi: storing the length of each event [-]
%   stormsi: storng the start ind of each storm event [-]
%
%
% inheritted from AWE-GEN-2D func findRadarEvents() with a bit changes to enable
% customizing seperation hour & radar resolution
%
% @ Yuting Chen (based on Peleg Code.)
% yuting.chen17@imperial.ac.uk
% Imperial College London
%
% Example: 
%    [lengthi,stormsi] = seperateRadarEvents( WAR, 0.02, 2, 60);
%    WAR(WAR < 0.02) = 0;
%    
%    hold on
%    bar((1:length(WAR))/12,WAR);
%    plot(stormsi/12,(stormsi>0)-1,'r.','markersize',50);
%    plot((stormsi+lengthi-1)/12,(stormsi>0)-1,'k.','markersize',50);
%    hold off
% -------------------------------------------------------------------------


arguments
    radarWAR double
    thresholdWAR double
    sepHour double = 2
    thresholdTime (1,1) double = 60
    radarReso (1,1) double = 5
end


%% Initilazing
WAR(radarWAR>=thresholdWAR,1)=1;
if radarReso == 5
    for i=6:size(WAR,1)-6 % Correcting time series of radar data for wet events hiatus
        if sum(WAR(i-5:i-1))>=1 && sum(WAR(i+1:i+5))>=1 && WAR(i)==0
            WAR(i)=1;
        end
    end
end

[lengthi,stormsi] = deal([]);

resoMin = 1/24/round((60/radarReso));
dryV = 0;
i = 1;

while(i<length(WAR))
    dur = 0;
    if WAR(i) <= dryV
        i = i+1;
    else
        % starts from i; to compute the duration of this storm
        dur = 0;
        starti = i;
        while i+dur<length(WAR) && WAR(i+dur)> dryV
            dur = dur + round(sepHour/(resoMin*24));% consectutive dry threshold is 2h;
        end
        % it has been dry days after sepHour h;
        if starti+dur<length(WAR)
            for j = dur:-1:1
                if WAR(i+j) > dryV 
                    dur = j+1;
                    break;
                end
            end
            i = i+dur;
        else
            disp([num2str(i),' ',num2str(length(WAR))]);
            i = length(WAR);
        end
        
        lengthi = [lengthi,i-starti];
        stormsi = [stormsi,starti];
    end
end


%% Checking for event length treshold

checkTime=lengthi<thresholdTime/radarReso;

for i=length(checkTime):-1:1
    if checkTime(i)==1
        stormsi(i) = [];
        lengthi(i) = [];
    end
end

end