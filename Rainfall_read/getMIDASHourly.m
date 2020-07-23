Year = 2006:2020;
% Arrange raw data and save it; get station information;
step = 2;
openfp = 'H:\DATA_CEDA\Hourly Rainfall Raw\';%'midas_rainhrly_',num2str(Year(fileInd)),'01-',num2str(Year(fileInd)),'12.txt']);
savefp = ['H:\DATA_CEDA\STATION_hourly_2006_2020\'];%['station',num2str(src_id(siteI)),'.txt'];
mkdir(savefp)
ArrangeRawData(openfp,savefp,Year,step,'hourly');
%%
A = load('H:\DATA_CEDA\STATION_hourly_2006_2020\station708.txt');
ob_start_time = datetime(datevec(A(:,1)));
ob_end_time = datetime(datevec(A(:,2)));
prcp_amt = A(:,3);

Heathrow = table(ob_start_time,ob_end_time,prcp_amt);
writetable(Heathrow,'H:\DATA_CEDA\STATION_hourly_2006_2020\station708.csv');


%%
function [] = ArrangeRawData(openfp,savefp,Year,step,resolution)
%% Data processing;
% import data of several files and transform the coordinate;
%% process raw date and save data of each station in .txt file;
switch(resolution)
    case 'hourly'
        for i = 1:step:length(Year)
            if i < length(Year)
                ArrangeHourlyData(openfp,savefp,Year(i:i+step-1));
                disp(['Saved! Year',num2str(Year(i)),'to',num2str(Year(i+step-1)),': Finished!']);
            else
                ArrangeHourlyData(openfp,savefp,Year(i));
                disp(['Saved! Year',num2str(Year(i)),': Finished!']);
            end
        end
    case 'subhourly'
        for i = 1:step:length(Year)
            if i < length(Year)
                ArrangeSubhourlyData(openfp,savefp,Year(i:i+step-1));
                disp(['Saved! Year',num2str(Year(i)),'to',num2str(Year(i+step-1)),': Finished!']);
            else
                ArrangeSubhourlyData(openfp,savefp,Year(i));
                disp(['Saved! Year',num2str(Year(i)),': Finished!']);
            end
        end
    otherwise
        error('Please Check!----ArrangeRawData(...resolution...)----@Yuting Chen');
end
end
function [] = ArrangeHourlyData(openfp,savefp,year)
% read rainfall data file;
%openfp = '\\icnas2.cc.ic.ac.uk\yc8316\Desktop\Data_CEDA\Hourlydata\';%'midas_rainhrly_',num2str(Year(fileInd)),'01-',num2str(Year(fileInd)),'12.txt']);
%savefp = ['D:','\DataInDifferentStation\'];%['station',num2str(src_id(siteI)),'.txt'];
%Year = [1980:1981];
N = length(year);
ob_end_time = [];
id_type = [];
ob_hour_count = [];
version_num = [];
src_id = [];
prcp_amt = [];
prcp_dur = [];
%id = [];
%met_domain_name = [];
%rec_st_ind = [];
%prcp_amt_q = [];
%prcp_dur_q = [];
%meto_stmp_time = [];
%midas_stmp_etime = [];
%prcp_amt_j = [];
for fileInd = 1:N
    fid = fopen([openfp,'midas_rainhrly_',num2str(year(fileInd)),'01-',num2str(year(fileInd)),'12.txt']);
    C = textscan(fid,'%s %d %s %d %d %s %d %d %f %d %f %d %s %s %s','delimiter',',','emptyvalue',NaN);
    fclose(fid);
    ob_end_time = [ob_end_time;C{1}];
    id_type = [id_type;C{3}];%'','CLBR','RAIN'
    ob_hour_count = [ob_hour_count;C{4}];%Observation hour count
    version_num = [version_num;C{5}];%Use the row with '1', as this has been quality checked by the Met Office
    src_id = [src_id;C{7}];
    prcp_amt = [prcp_amt;C{9}];%Precipitation amount, Units = 1mm, reported to the nearest 0.1 mm
    prcp_dur = [prcp_dur;C{10}];%Precipitation duration (<24 hr) minutes
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% note: in some data, there is no useful info from prcp_dur.
    %id = [id;C{2}];%Raingauge number
    %met_domain_name = [met_domain_name;C{6}];
    %rec_st_ind = [rec_st_ind;C{8}];
    %prcp_amt_q = [prcp_amt_q;C{11}];
    %prcp_dur_q = [prcp_dur_q;C{12}];
    %meto_stmp_time = [meto_stmp_time;C{13}];
    %midas_stmp_etime = [midas_stmp_etime;C{14}];
    %prcp_amt_j = [prcp_amt_j;C{15}];
    clear fid;
    
end
clear C

% filter data with low quality;
lowQ_index = find(version_num == 0);
ob_end_time(lowQ_index) = [];
id_type(lowQ_index) = [];
ob_hour_count(lowQ_index) = [];
version_num(lowQ_index) = [];
src_id(lowQ_index) = [];
prcp_amt(lowQ_index) = [];
prcp_dur(lowQ_index) = [];

% filter data with long observation duration. (> 2h)
longT_index = find(ob_hour_count > 2);
ob_end_time(longT_index) = [];
id_type(longT_index) = [];
ob_hour_count(longT_index) = [];
version_num(longT_index) = [];
src_id(longT_index) = [];
prcp_amt(longT_index) = [];
prcp_dur(longT_index) = [];


% read src_id file and locate the station;
%{
ILatLon = xlsread('src_ids.xlsx', 'CYT','A2:C17339');
[B,I] = sort(ILatLon(:,1));
ILatLon = ILatLon(I,:);
stationId = double(unique(src_id));
id_lat_lon = [stationId,zeros(numel(stationId),2)];
%find corresponding location of each site;
for i = 1:numel(stationId)
    ii = find(ILatLon(:,1)==stationId(i));
    if ~isempty(ii)
        id_lat_lon(i,2:3) = ILatLon(ii,2:3);
    end
end
%}

% ILatLon = xlsread('station_detailsForHourlyData1995&2017.xlsx','excel_list_station_details','A2:C352');
stationId = double(unique(src_id));
id_latLon_EN = [stationId,zeros(numel(stationId),4)];

% plot(id_latLon_EN(:,4),id_latLon_EN(:,5),'.');
ob_end_time = datenum(datetime(ob_end_time, 'InputFormat', 'yyyy-MM-dd HH:mm'));%1day-->1
ob_start_time = ob_end_time - double(ob_hour_count)/24;
% ts_intensity = nan(length(ob_hour_count),1);% mm/h

% calculate the intensity (depth in 1 hour);
% normal_index = find(ob_hour_count == 1);
% ts_intensity(normal_index) = prcp_amt(normal_index);
% special_index = find(ob_hour_count > 1);
% ts_intensity(special_index) = prcp_amt(special_index)./double(ob_hour_count(special_index));

%% save time series of each station in different file;

    index = find(src_id == 708);
    fileName = ['station',num2str(708),'.txt'];
    v1 = [ob_start_time(index),ob_end_time(index),prcp_amt(index)];%ts_intensity(index)];
    % v2 = id_type(index);
    save([savefp,fileName],'v1','-ascii','-double','-append');


% for siteI = 1:length(id_latLon_EN)
%     % save ob_start_time ob_end_time ts_intensity id_type for each src_id;
%     index = find(src_id == id_latLon_EN(siteI,1));
%     fileName = ['station',num2str(id_latLon_EN(siteI,1)),'.txt'];
%     v1 = [ob_start_time(index),ob_end_time(index),prcp_amt(index)];%ts_intensity(index)];
%     % v2 = id_type(index);
%     save([savefp,fileName],'v1','-ascii','-append');
% end
end

function [] = ArrangeSubhourlyData(openfp,savefp,year)
% read rainfall data file;
N = length(year);
ob_end_time = [];
src_id = [];
prcp_amt = [];
for fileInd = 1:N
    fid = fopen([openfp,'midas_rainsub_',num2str(year(fileInd)),'01-',num2str(year(fileInd)),'12.txt']);
    C = textscan(fid,'%s %s %f %s %f %f %f %s %f','delimiter',',','emptyvalue',NaN);
    fclose(fid);
    ob_end_time = [ob_end_time;C{1}];
    src_id = [src_id;C{5}];
    prcp_amt = [prcp_amt;C{6}];%Precipitation amount, Units = 1mm, reported to the nearest 0.001 mm?
    
    clear fid;
    
end
clear C

% ILatLon = xlsread('station_detailsForHourlyData1995&2017.xlsx','excel_list_station_details','A2:C352');
stationId = double(unique(src_id));
id_latLon_EN = [stationId,zeros(numel(stationId),4)];

% plot(id_latLon_EN(:,4),id_latLon_EN(:,5),'.');
ob_end_time = datenum(datetime(ob_end_time, 'InputFormat', 'yyyy-MM-dd HH:mm'));%1day-->1

ts_depth = prcp_amt;


%% save time series of each station in different file;
for siteI = 1:length(id_latLon_EN)
    % save ob_start_time ob_end_time ts_intensity id_type for each src_id;
    index = find(src_id == id_latLon_EN(siteI,1));
    fileName = ['station',num2str(id_latLon_EN(siteI,1)),'.txt'];
    v1 = [ob_end_time(index),ts_depth(index)];
    %v2 = id_type(index);
    save([savefp,fileName],'v1','-ascii','-append');
end
end