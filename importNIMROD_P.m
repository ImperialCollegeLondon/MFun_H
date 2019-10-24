function [PRS,status] = importNIMROD_P(XX,YY,YEARRANGE,varargin)
% IMPOERTNIMROD_P return the time series of the precipitation over the
% space
%
% Radar Format: NIMROD weather radar data in UK.
% Default radar file location: 'D:\UK_Radar_Matlab'
%
%
% Input: XX: single value/vector/matrix <double>
%        YY: same size as XX.
%        YEARRANGE:
% Output:PRS: 3D Matrix<int16> scale-32;

% Example 1:
%     %Import PRS data at Yare station
%     XX = 618200;
%     YY = 308200;
%     YEAR = 2007:2008;
%     [PRS,status] = importNIMROD_P(pointX,pointY,YEAR)
%     status


% Example 2:
%     % Import PRS data for one catchment
%     CatchName = 'Lamarsh';
%     sfName = [CatchName,'_TOPO_DATA_500CellSize.mat'];
%     C = load(sfName);
%     [XX,YY] = meshgrid(C.y_yr,C.x_yr);
%     YEAR = 2007:2007;
%     [PRS,status] = importNIMROD_P(XX,YY,YEAR);
%     status

% by Yuting CHEN
% Imperial College London
% yuting.chen17@imperial.ac.uk


status = 0;
rp = 'E:\UK_Radar_2007_2008\';
commandwindow
fprintf('Now the file path is: %s\n',rp);

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
startD = datetime(min(YEARRANGE),1,1);
endD = datetime(max(YEARRANGE),12,31);

%%
tic

totalDays = datenum(endD)-datenum(startD);

nan_day = 0;


PRS = NaN(size(XX,1),size(XX,2),24*12*(totalDays+1),'single');

for ind_day = 0:totalDays
    tic
    da = startD + ind_day;
    filename = [num2str(da.Year),'_',num2str(da.Month),'_',num2str(da.Day),'.mat'];
    try
        A = load([rp,filename]);
    catch
        try
            rp = 'D:\UK_Radar_Matlab\';
            A = load([rp,filename]);
        catch
            try
            rp = 'H:\UK_Radar_Matlab\';
            A = load([rp,filename]);
            catch
                fprintf('Check!');
            end
        end
        
    end
    DAT = A.DAT;
    times = fieldnames(DAT,'-full');

    for minu = 1:length(times) % for each day; resolution: 5min;
        eval(['aaa = DAT.' times{minu} ';']);
        if ~isnan(aaa.rr)         
            dx = aaa.rl_gen_hd(6);
            dy = aaa.rl_gen_hd(4);
            x0 = aaa.rl_gen_hd(5);
            y0 = aaa.rl_gen_hd(3);
            loci = @(Y)round((y0-Y)/dy+1);
            locj = @(X)round((X-x0)/dx+1);
            sei = loci(YY);
            sej = locj(XX);
            RR = aaa.rr;
            if sum(size(RR)==[2175,1725])==2
                rain = reshape(RR(unique(sei),unique(sej)),size(sei));% same matrix as in XX and YY;
                PRS(:,:,ind_day*24*12 + minu) = rain;
            else
                rain = NaN(size(XX,1),size(XX,2),'single');
                PRS(:,:,ind_day*24*12 + minu) = rain;
            end
        else
            % save it as NaN number when the data is not available
            rain = NaN(size(XX,1),size(XX,2),'single');
            PRS(:,:,ind_day*24*12 + minu) = rain;
            nan_day = nan_day + 1;
            disp([times{minu},' nan+1-->',num2str(nan_day)]);
        end
    end
    disp(['Iter: Day:',num2str(ind_day)]);
    toc
%     image(double(aaa.rr));hold on;
%     %contourf(sei,sej,RR(unique(sei),unique(sej)),'r.');hold on;
%     plot(sej(:),sei(:),'r.');
%     pause(0.05);
    
end
toc

disp('finished');

PRS = int16(PRS);

end

%%%%%%%%%%%%%%%%%%%%%%
% See output quality;
% for i = 1:size(PRS,3)
%     imagesc(squeeze(PRS(:,:,i)));
%     pause(0.02)
% end



