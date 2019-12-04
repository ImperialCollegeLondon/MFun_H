


function [t_storm,p_storm,st_storm] = deriveStormStructure(P,dat)

%% compute the statistical properties ralated to storm.
%
% Input:   P:cell(stationNum)
%          dat: cell(stationNum) datenum
%
% Output:  t_storm:cell(1,stationNum), save the duration of each storm;
%          p_storm:cell(1,stationNum), save the total depth of each storm;
%          st_storm:cell(1,stationNum), save the start time of each storm;
%
% Example: simSeriesPath = ['1_hour_Simdata_1to11_EA.mat'];
%          dt = 1;
%          [SIM_P,SIM_dat] = importSimulatedData(simSeriesPath,dt,IND);
%          [t_storm,p_storm,st_storm] = deriveStormStructure(SIM_P,SIM_dat);
%          sim_storm = struct('t_storm',t_storm,'p_storm',p_storm,'st_storm',st_storm);
%
% Author Yuting Chen
%
stationNum = length(P);
t_storm = cell(1,stationNum);
p_storm = cell(1,stationNum);
st_storm = cell(1,stationNum);

reso = dat{1}(2)-dat{1}(1);
dryV = 0.2;%1/24;
for sta = 1:stationNum
    i = 1;
    while(i<length(P{sta}))
        dur = 0;
        if P{sta}(i) < dryV
            i = i+1;
        else
            % starts from i; to compute the duration of this storm
            dur = 0;
            starti = i;
            while i+dur<length(P{sta}) && P{sta}(i+dur)>= dryV
                dur = dur + round(12/(reso*24));% consectutive dry threshold is 6h;
            end
            % it has been dry days after 6h;
            if starti+dur<length(P{sta})
                for j = 1:dur
                    if P{sta}(i+j) < 1/24%dryV
                        dur = j;
                        break;
                    end
                end
                i = i+dur;
            else
                disp([num2str(i),' ',num2str(length(P{sta}))]);
                i = length(P{sta});
            end
            
            t_storm{sta} = [t_storm{sta},dat{sta}(i)-dat{sta}(starti)];
            p_storm{sta} = [p_storm{sta},sum(P{sta}(starti:i))];
            st_storm{sta} = [st_storm{sta},dat{sta}(starti)];
        end
    end
end
end
