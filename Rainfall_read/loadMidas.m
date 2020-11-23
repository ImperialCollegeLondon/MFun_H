function [T,P,X,Y,dt] = loadMidas(Region,DFolder)
% LOADMIDAS will extract the MIDAS hourly rainfall data for a specific
% region;
%
% ---------------------------------------------- %
% Might Need to be rewritten into your data
% EXAMPLE:
% A = load('CompleteTS.mat');
% T = A.CompleteTS.T{1};
% P = A.CompleteTS.P;
% X = A.CompleteTS.ENA(:,1);
% Y = A.CompleteTS.ENA(:,2);
% dt = 1;%unique(round(diff(T)*24*60))/60;
% Might Need to be rewritten into your data %

sfp = 'H:\Data_CEDA\DataInDifferentStation_New discardLongObservationCount\';
stationInfo_fp = 'E:\OneDrive - Imperial College London\MSc_Dissertation\STOCHASTIC_MODEL0528\StationInfo.mat';
dreamYear = 2005:2017;%1990:2017;

    
load(stationInfo_fp,'StationInfo');

[Xrange,Yrange] = getXYrange(Region);

[src_id] = getSrcId(Xrange,Yrange,StationInfo,dreamYear);

[T,P,X,Y,dt] = readTS(src_id,StationInfo,sfp,dreamYear);

[T,P,X,Y] = formatOutput(T,P,X,Y);

save(sprintf('%sTS_obs.mat',DFolder),'T','P','X','Y','dt','src_id');

end

function [Xrange,Yrange] = getXYrange(Region)
switch(Region)
    case 'Birmingham'
        % Midlands
        Xrange = [305000,605000]/1000;% 200km
        Yrange = [187000,387000]/1000;% 200km
    case 'Plynlimon'
        % Plynlimon
        Xrange = [235000,331000]/1000;
        Yrange = [240000,330000]/1000;
    case 'Stour'
        % Stour
        Xrange = [555000,610000]/1000;
        Yrange = [225000,265000]/1000;
        
end
end

function [src_id] = getSrcId(Xrange,Yrange,StationInfo,dreamYear)

Year = 1957:2017;
% dreamYear = [1990,2017];

src_id = [];% store src_id;

for i = 1:length(StationInfo.src_id)
    if isempty(find(StationInfo.yearMatrix(i,dreamYear(1)-Year(1)+1:dreamYear(end)-Year(1)+1)==0, 1))
        if StationInfo.E(i)<Xrange(2) && StationInfo.E(i)>Xrange(1)&& ...
                StationInfo.N(i)>Yrange(1) && StationInfo.N(i)<Yrange(2)
            
            src_id = [src_id,StationInfo.src_id(i)];
            
        end
    end
end


end

function [T,P,X,Y,dt] = readTS(src_id,StationInfo,sfp,dreamYear)

dt = 1;

[X,Y] = getCoor(src_id,StationInfo);
[T, P] = create_gauge_series_hourly( sfp,src_id,dreamYear);

end

function [X,Y] = getCoor(src_id,StationInfo)

id = arrayfun(@(id)find(StationInfo.src_id == id),src_id);
X = StationInfo.E(id);
Y = StationInfo.N(id);

end

function [ gauge_dat, gauge_series, station_id] = create_gauge_series_hourly( rp,station_id,time_domain)
% Input Parameters:
%       rp: root path;
%       station_id: vertex of station_id;
%       time_domain: vertex,[minYear, maxYear];

% Note:
%       Data should be continuous with the time_domain;
%       guage_dat is vertex;
%       gauge_series is cell;
%       generated gauge_dat would be *continuous* for hourly scale;
%% Load data
% search_domain_x = [681.193-64 681.193+64];
% search_domain_y = [237.593-64 237.593+64]; % around albis


ind = numel(station_id);

depth = cell(1,ind);
date_raw = cell(1,ind);

for i = 1:ind
    fid = fopen([rp,'station',num2str(station_id(i)),'.txt']);
    reso = 60;
    C = textscan(fid,'%f %f %f','emptyvalue',NaN);
    dat = round(24*C{2})/24; 
    p = C{3};
    fclose(fid);
    
    depth{i}=p(~isnan(p));
    date_raw{i} = dat(~isnan(p));
    
    if ~isnan(time_domain)
        index = find(date_raw{i}<datenum(time_domain(end)+1,1,1) & date_raw{i}>=datenum(time_domain(1),1,1));
        date_raw{i} = date_raw{i}(index);
        depth{i} = depth{i}(index);
    end
end
for i = 1:ind
    
    [temp_date,IA,~] = unique(date_raw{i});

    [depth{i},date_raw{i}] = interpRain(temp_date,depth{i}(IA),maxMin(date_raw),minMax(date_raw),reso);
    
end
clear C dat;

gauge_dat = date_raw{1};%datevec(date_num(:));
gauge_series = depth;

    function res = minMax(datacell)
        % to find the common min value in datacell{i}
        temp = [];
        for i1 = 1:length(datacell)
            temp = [temp,max(datacell{i1})];
        end
        res = min(temp);
    end

    function res = maxMin(datacell)
        temp = [];
        for i2 = 1:length(datacell)
            temp = [temp,min(datacell{i2})];
        end
        res = max(temp);
        
    end
    function [depth,date] = interpRain(oridate,oridepth,minDate,maxDate,reso)
        
        maxdiff = max(diff(oridate));
        if maxdiff>10
            ind2 = find(diff(oridate) == maxdiff);
            datetime(datevec(oridate(ind2))) %#ok<FNDSB>
            sprintf('Pay Attention! poor quality data: max diff:%dday',maxdiff);
        end
        
        %create series for calibration process
        depth = interp1(oridate,oridepth,[minDate:1/24/(60/reso):maxDate]');
        date = [minDate:1/24/(60/reso):maxDate]';
        threshold = 0;%0.2;% daily threshold;
        depth(depth<=threshold/24/(60/reso)) = 0;
        
    end
end

function [T,P,X,Y] = formatOutput(T,P,X,Y)
T = reshape(T,1,[]);
P = cellfun(@(x)reshape(x,1,[]),P,'UniformOutput', false);
P = reshape(P,1,[]);
X = reshape(X,[],1);
Y = reshape(Y,[],1);
end
