function [status] = aggregate_NIMROD_Hour(YEAR,varargin)
% AGGREGATENIMROD2HOUR return the time series of the precipitation over the
% space
%
% Here we are processing 5min NIMROD Radar Composite which we got from UK
% Met Office (Thanks for supper powerful BADC!)
%
% Radar Format: NIMROD weather radar data in UK.
% Default radar file location: 'D:\UK_Radar_Matlab'
%
%
% Input: YEARRANGE:
% Output:PRS: 3D Matrix<int16> scale-32;
%
%
% Example 1:
%     %Import PRS data at Yare station
%     YEAR = 2007;
%     [PRS,status] = aggregateNIMROD2Hour(YEARRANGE)
%     status
%
% by Yuting CHEN
% Imperial College London
% yuting.chen17@imperial.ac.uk
% Update: 2019.12.02


status = 0;
commandwindow
rp = 'H:\DATA_RADAR\UK_Radar_Matlab\';
fprintf('Now the file path is: %s\n',rp);


if ~isscalar(YEAR)
    error('Check Input:YEAR');
end


% s = input('Make sure it is right? Y/N [Y]:','s');
s = 'Y';

if strcmp(s,'Y')
    %%% okay;
elseif strcmp(s,'N')
    fprintf('Please Check Path\n');
    status = 1;
    return
else
    fprintf('Please Check Path\n');
    status = 1;
    return
end

% yearRange = 2008:2009;
startD = datetime(YEAR,1,1);
endD = datetime(YEAR,12,31);

%%
tic

totalDays = datenum(endD)-datenum(startD);

% PRS = NaN(2175,1725,12*(totalDays+1),'single');
fid = fopen(['Year',num2str(YEAR),'_dialogFile.txt'],'w');

for ind_day = 0:totalDays
    tic
    da = startD + ind_day;
    filename = [num2str(da.Year),'_',num2str(da.Month),'_',num2str(da.Day),'.mat'];
    try
        load([rp,filename],'DAT');
    catch
        continue
    end
    times = fieldnames(DAT,'-full');
    c1 = 2175;
    c2 = 1725;
    PRS0 = NaN(c1*c2,24*60/5);%2D data (x*y,time)
    
    for minu = 1:length(times) % for each day; resolution: 5min;
        aaa = getfield(DAT,times{minu});
        if ~isnan(aaa.rr)         
            dx = aaa.rl_gen_hd(6);
            dy = aaa.rl_gen_hd(4);
            x0 = aaa.rl_gen_hd(5);
            y0 = aaa.rl_gen_hd(3);
            if sum(size(aaa.rr)==[c1,c2])==2
                PRS0(:,minu) = reshape(aaa.rr,[],1);
            else
                fprintf(fid,'%s\n',times{minu});
            end
        else
            % save it as NaN number when the data is not available
            fprintf(fid,'%s\n',times{minu});  
        end
    end
    disp(['Iter: Day:',num2str(ind_day)]);
    clear DAT
    
    % aggregate to hourly scale.
    PRS0 = int16(imresize(PRS0,'Scale',[1,1/12]));%[c1 c2 24]));
    save(['H:\DATA_RADAR\UK_Radar_HourlyAggregate\',filename],'PRS0');
    PRS0 = 0;
    toc
    % image(squeeze(aa(:,:,1)));hold on;
    % pause(0.05);
    
end
toc

disp('finished');


end

%%%%%%%%%%%%%%%%%%%%%%
% See output quality;
% for i = 1:size(PRS,3)
%     imagesc(squeeze(PRS(:,:,i)));
%     pause(0.02)
% end



