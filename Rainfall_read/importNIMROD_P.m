function  [DATA,status] = importNIMROD_P(XX,YY,TIME,varargin)
% 
% IMPOERTNIMROD_P return the time series of the precipitation over the
% space
%
% Radar Format: NIMROD weather radar data in UK.
% Default radar file location: 'D:\UK_Radar_Matlab'
%
%
% Input: XX: single value/vector/matrix <double>
%        YY: same size as XX.
%        TIME: YEARRANGE <double> or DATATIME<datetime>
% Output:PRS: 3D Matrix<int16> scale-32;
%
% Example 1:
%     %Import PRS data at Yare station
%     XX = 618200;% Easting
%     YY = 308200;% Northing
%     YEAR = 2007:2008;
%     [PRS,status] = importNIMROD_P(pointX,pointY,YEAR)
%     status
%
%
% Example 2:
%     % Import PRS data for one catchment
%     x_yr = 1000:1000:5000;
%     y_yr = 10000:-1000:2000;
%     [XX,YY] = meshgrid(x_yr,y_yr);
%     YEAR = 2007:2007;
%     [PRS,status] = importNIMROD_P(XX,YY,YEAR);
%     status
%
% by Yuting CHEN
% Imperial College London
% yuting.chen17@imperial.ac.uk
%
% Note: OUTPUT format: [loc1,loc2,time]
% Exactly in the same acsending/descending/1stNorthing/1stEasting order as input XX, YY.
% 
% Update:
%        2020.05.22 Change output into <class>
% 

status = 0;
rp = 'K:\UK_Radar_Matlab\';
fprintf('Now the file path is: %s\n',rp);

if isdatetime(TIME)
    DATATIME = TIME;
    startTime = datetime(datevec(DATATIME(1)));
    startD = datetime(startTime.Year,startTime.Month,startTime.Day);
    endTime = datetime(datevec(DATATIME(end)));
    endD = datetime(endTime.Year,endTime.Month,endTime.Day);
else
    YEARRANGE = TIME;
    startD = datetime(min(YEARRANGE),1,1);
    endD = datetime(max(YEARRANGE),12,31);
end
PTime = startD:minutes(5):endD+days(1)-minutes(5);
%%
tic

totalDays = datenum(endD)-datenum(startD);

nan_day = 0;
% [sei,sej] = deal(NaN);

PRS = -1*ones(size(XX,1),size(XX,2),24*12*(totalDays+1),'int16');

for ind_day = 0:totalDays
    [sei,sej] = deal(NaN);
    tic
    da = startD + ind_day;
    filename = [num2str(da.Year),'_',num2str(da.Month),'_',num2str(da.Day),'.mat'];
    da.Format = 'yyyyMMddHHmm';
    try
        load([rp,filename],'DAT');
        
        times = fieldnames(DAT,'-full');
        
        rains = struct2cell(structfun(@(R)getRain(R), DAT, ...
            'UniformOutput', false));
        for minu = 1:length(times)
            PRS(:,:,ind_day*24*12 + minu) = rains{minu};
        end
        
    catch me
        me
    end
    disp(['Iter: Day:\n',num2str(ind_day)]);
    
end
toc

disp('finished');

PRS = int16(PRS);

DATA = RainfallDataClass(PRS,-1,32,'mm/h',PTime',XX,YY,'');

if isdatetime(TIME)
    [isInPeriod,~] = ismember(PTime',TIME);
    DATA = extractOnePeriod(DATA,isInPeriod);
end

%%
    function rain = getRain(aaa)
        if isnan(sei)
            dx = aaa.rl_gen_hd(6);
            dy = aaa.rl_gen_hd(4);
            x0 = aaa.rl_gen_hd(5);
            y0 = aaa.rl_gen_hd(3);
            loci = @(Y)round((double(y0)-double(Y))/double(dy)+1);
            locj = @(X)round((double(X)-double(x0))/double(dx)+1);
            sei = loci(YY);
            sej = locj(XX);
        end
        
        if ~isnan(aaa.rr)
            RR = aaa.rr;
            if sum(size(RR)==[2175,1725])==2
                rain = RR(sub2ind(size(RR),sei,sej));
            elseif sum(size(RR) == [775,640])==2
                rain = RR(sub2ind(size(RR),sei,sej));
            else
                rain = ones(size(XX,1),size(XX,2),'int16')*(-1);
            end
        else
            % save it as '-1' when the data is not available
            rain = ones(size(XX,1),size(XX,2),'int16')*(-1);
            nan_day = nan_day + 1;
        end
        
    end
end


