function [Prs_coarse,status,scale] = aggregateRADAR(Prs,dt,dx,T,ijSize,Thres_abn);
% AGGREGARERADAR finish the aggregation for radar data in space and in
% time.

% Note: In the case that T<dt or X<dx:
%         interpolation will not be excuted.
%         rectangular will be processed.

% Input:
%       Prs_RADAR: Radar time series in space; 3D Matrx<int16> X.*Y.*T
%       dt: original resolution in time; <double>
%       dx: original resolution in space; <double>
%       T: aggregated time resolution of interest; <double>
%       ijSize: spatial size after aggregation; vector<double>
%       Thres_abn: threshold for discarding the abnormal value;
%
% Output:
%       Prs_coarse: aggragetd Radar time series in space; 3D Matrix<int16> X.*Y.*T
%       status: diagonal status <0/1>
%       scale: scaling factor in Prs_coars; true value should be:
%              Prs_coarse/scale
% Example: 
%       radar = load(['PRS_radar_',num2str(year_radar),'_',CatchName,'.mat']);
%       topo = load(sfName);
%       dt = 5; % [min];
%       dx = 1; % [km];
%       T = 60; % [hour/min];
%       ijSize = [154 230];
%       Thres_abn = 200*32;% 200mm/h * 32;
%       
%       [Prs_coarse,status,Prs_scale] = aggregateRADAR(radar.PRS,dt,dx,T,ijSize,Thres_abn);
%
% @ Yuting CHEN

%%% Check Data;
if ~(dt <= T)
    fprintf(['Check Input Variables..\n',...
        'dt = %f T = %f\n',...
        'Normally: dt<T. \n',...
        'Are you sure it is correct?\n'],dt,dx);
    while(1)
        reply = input('Are you sure it is correct? Y/N','s');
        if strcmp(reply,'Y') || strcmp(reply,'y')
            break;
        elseif strcmp(reply,'N')|| strcmp(reply,'n')
            fprintf('Please Reinput\n');
            return;
        else
            %%%
        end
    end
end

%%% single/int16/double
Data = Prs;

%%% Processing For NaN Value.
if sum(isnan(Data(:))) ~= 0
    
    fprintf('*NaN value* Check input Radar data.\n');
    return;

end
%%% Processing For AbnormalBig Value.
Data(Data>Thres_abn) = NaN;
Data(Data<0) = NaN;

%%% Aggregate In Space / In Time.
% TRY in time:
Ti = size(Data,3);
try
    % in time: aggregate into 1 hour
    scale = round(T/dt);
    if scale ~= 1
        Prs_2 = NaN(size(Data,1),size(Data,2),ceil(Ti/round(T/dt)),'single');
        for i = 1:size(Data,1)
            for j = 1:size(Data,2)
                ts_aux = squeeze(Data(i,j,:));
                Prs_2(i,j,:) = aggregateTime(ts_aux,round(T/dt),'mean');
            end
        end
    else
        Prs_2 = Data;
    end
    
catch
    fprintf('Error happens in # aggregateRADAR.m \n #');
end


Data = Prs_2;
Ti = size(Data,3);

% First in space:
if size(ijSize,1)~=size(Data,1) || size(ijSize,2)~=size(Data,2)
    Prs_1 = NaN(ijSize(1),ijSize(2),Ti,'single');
    
    for i = 1:Ti
        
        dat_aux = squeeze(Data(:,:,i));
        mat_aux = imresize(dat_aux,ijSize,'nearest');% replicate without interpolation;
        
        Prs_1(:,:,i) = single(mat_aux);
        % --------------------------------------------------------------------%
        % Notice: Here, when doing disaggregation, some location discrepancies
        % might exist due to the different resolution but can be accpetable
        % (<500m).
        % Futher improvement: process the inner part and the
        % boundary vector seperately;
        % by. Yuting CHEN
        % --------------------------------------------------------------------%
    end
else
    Prs_1 = Data;
end


status = 1;


%%% Return Status. Prs_coarse.

Prs_coarse = Prs_1;


end






%% Auxillary function
function out=aggregateTime(data,scale,method)

A=buffer(data,scale);

if scale==1
    out=data;
else
    
    
    if (strcmp(method,'mean'))
        %             out=mean(A);
        out=nanmean(A);
    elseif(strcmp(method,'sum'))
        %             out=sum(A);
        out=nansum(A);
    else
        error('function out=aggregate() BUG 1 !!!!')
    end
    
end
% disp(length(out));
out=out(:);
end

