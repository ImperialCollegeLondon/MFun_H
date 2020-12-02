% extractRadar_demo.m
%
% this .m file shows how to extract radar info for a study area,
% from the 'F:\UK_Radar_Matlab\'.
%
% Yuting Chen
% Imperial College London

%% main
clear;clc

[PRS,T,XX,YY] = extractDays();

% plt time series at one location
figure;
plot(T,squeeze(PRS(1,1,:)));
xlabel('Time');
ylabel('Intensity [mm/h]');

%% auxillary func
function [PRS,T,XX,YY] = extractDays()
% Output:
%       PRS: intensity[mm/h] <double>
%       T: time <datetime>
%       XX: easting[m] <double>
%       YY: northing[m] <double>

%% configuration
filepath = ['F:\UK_Radar_Matlab\'];
x_yr = 530e3:1e3:540e3;
y_yr = 180e3:1e3:185e3;
[XX,YY] = meshgrid(x_yr,y_yr);
T_days = datetime(2010,8,21:22);

%% extraction

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

%% transformation
PRS = double(PRS);
PRS(PRS == -1) = NaN;
PRS = PRS/32;% unit: [mm/h]

%% nested func
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
