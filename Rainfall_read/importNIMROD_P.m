function  [DATA,status] = importNIMROD_P(XX,YY,TIME,varargin)
% IMPOERTNIMROD_P returns the time series of the precipitation over the
% space
%
% Radar Format: NIMROD weather radar data in UK.
% Default radar file location: 'K:\UK_Radar_Matlab'
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
%     YEAR = datetime(2000,1,1):datetime(2001:12,31,10,20,0);
%     [PRS,status] = importNIMROD_P(pointX,pointY,YEAR)
%     status
%
%
% Example 2:
%     % Import PRS data for one catchment
%     x_yr = 530e3:1e3:540e3;
%     y_yr = 180e3:1e3:185e3;
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
    t_temp = datetime(datevec(DATATIME(1)));
    startD = datetime(t_temp.Year,t_temp.Month,t_temp.Day);
    t_temp = datetime(datevec(DATATIME(end)));
    endD = datetime(t_temp.Year,t_temp.Month,t_temp.Day);
else
    YEARRANGE = TIME;
    startD = datetime(min(YEARRANGE),1,1);
    endD = datetime(max(YEARRANGE),12,31);
end
T_days = startD:days(1):endD;
PrsTime = startD:minutes(5):endD+days(1)-minutes(5);

%%
tic

T = datetime(T_days(1)):minutes(5):datetime(T_days(end)+days(1)-minutes(5));
PRS = int16(-ones([size(XX),numel(T)]));
itag = 0;

for oneDay = T_days
    
    filename = [num2str(oneDay.Year),'_',num2str(oneDay.Month),'_',num2str(oneDay.Day),'.mat'];
    oneDay.Format = 'yyyyMMddHHmm';
    
    [loci,locj] = deal([]); % relocate for each day
    
    try
        load([filepath,filename],'DAT');
        times = fieldnames(DAT,'-full');
        
        rains = struct2cell(structfun(@(r)getRain(r,XX,YY), DAT, ...
            'UniformOutput', false));
        rains = permute(cell2mat(rains),[2,3,1]);
        PRS(:,:,days(1)/minutes(5)*itag+(1:length(times))) = rains;
        fprintf('%s finished\n',setfield(oneDay,'Format','yyyy-MM-dd'));
    catch
        % no file
        fprintf('%s no data\n',setfield(oneDay,'Format','yyyy-MM-dd'));
        PRS(:,:,days(1)/minutes(5)*itag+(1:length(times))) = -1;
    end
    
    itag = itag+1;
    
end
toc

disp('finished');

PRS = int16(PRS);

DATA = RainfallDataClass(PRS,-1,32,'mm/h',PrsTime',XX,YY,'');

if isdatetime(TIME)
    [isInPeriod,~] = ismember(PrsTime',TIME);
    DATA = extractOnePeriod(DATA,isInPeriod);
end

%%
    function rain = getRain(oneSnapshot,XX,YY)
        
        if isempty(loci) % prevent repeated locating during a same day
            dx = oneSnapshot.rl_gen_hd(6);% 1000
            dy = oneSnapshot.rl_gen_hd(4);% 1000
            x_ul = oneSnapshot.rl_gen_hd(5);% ul: upper left
            y_ul = oneSnapshot.rl_gen_hd(3);
            func_loci = @(Y)round((double(y_ul)-Y)/double(dy)+1);
            func_locj = @(X)round((X-double(x_ul))/double(dx)+1);
            loci = func_loci(YY);
            locj = func_locj(XX);
        end
        
        if ~isnan(oneSnapshot.rr)
            RR = oneSnapshot.rr;
            if sum(size(RR)==[2175,1725])==2 | sum(size(RR) == [775,640])==2 %#ok<OR2>
                rain = RR(sub2ind(size(RR),loci,locj));
            else
                rain = ones(size(XX,1),size(XX,2),'int16')*(-1);
            end
        else
            % no data
            rain = ones(size(XX,1),size(XX,2),'int16')*(-1);
        end
        
        rain = reshape(rain,[1,size(rain)]);
        
    end
end


