function combineSim_batch(fileNameSeries,year,SCALE,goalFileName)
% COMBINESIM_BATCH combine a series of sim files (runoff) into one file
% INPUT
%     fileNameSeries
%     year
%     goalFileName
% OUTPUT
%     ETR_val
%     gwDisPart_val
%     leakage_val
%     O
%     Qpoint_val
%     QpointS_val
%     Runoff_time
%     Runoff_val
%     Vf_val
% EXAMPLE
% year = 2007:2015;
% fileNameSeries = [];
% expNo = 5;
% CatchName = 'Welland';
% ii = 1;
% for YEAR = year
%     fileNameSeries{ii} = ['sim_batch\Radar_GEAR_',CatchName,'\output_',CatchName,...
%         '_8736for_200m_9Years',num2str(expNo),'_',...
%         sprintf('%4d',YEAR),'.mat'];
%     try
%         load(fileNameSeries{ii});
%     catch
%         fileNameSeries{ii} = ['sim_batch\Radar_GEAR_',CatchName,'\output_',CatchName,...
%             '_8760for_200m_9Years',num2str(expNo),'_',...
%             sprintf('%4d',YEAR),'.mat'];
%     end
%     ii = ii+1;
% end
% SCALE = 1;
% goalFileName = ['F:\data_sim_NS\',CatchName,'_9Years',num2str(expNo),'_runoff.mat'];
% combineSim_batch(fileNameSeries,year,SCALE,goalFileName)
%
% NOTE
% 
% UPDATE
%     19 Dec 2019.


% NORMALIZATION
if ~iscell(fileNameSeries)
    fileNameSeries = {fileNameSeries};
end

% INITIALIZTION
fileNum = length(fileNameSeries);
dt = 1; % hour
h2sec = 3600;

Runoff_time = datenum(year(1),1,1):SCALE/24:datenum(year(end),12,31,23,0,0); %#ok<NASGU>
[Runoff_val,Qpoint_val,QpointS_val,Vf_val,gwDisPart_val,ETr_val,leakage_val,O] = deal([]);

% COMBINE
itag = 1;
for i = 1:fileNum
    % LOAD FILE
    load(fileNameSeries{i},'OUTPUT','TOPO_DAT');
    
    % ----- LOCATE RUNOFF
    for pointi = 1:size(OUTPUT.QpointC,2)
        
        ZZ = OUTPUT.QpointC(2:end,pointi)/1000*(TOPO_DAT.cellsize^2)/h2sec;
        ZZ = addVector(ZZ,year(i),SCALE,dt);
        Runoff_val{pointi}(itag:itag+numel(ZZ)-1) = ZZ;
        
        ZZ = OUTPUT.Qpoint(2:end,pointi)/1000*(TOPO_DAT.cellsize^2)/h2sec;
        ZZ = addVector(ZZ,year(i),SCALE,dt);
        Qpoint_val{pointi}(itag:itag+numel(ZZ)-1) = ZZ;
        
        ZZ = OUTPUT.QpointS(2:end,pointi)/1000*(TOPO_DAT.cellsize^2)/h2sec;
        ZZ = addVector(ZZ,year(i),SCALE,dt);
        QpointS_val{pointi}(itag:itag+numel(ZZ)-1) = ZZ;
        
    end
    % ----- LOCATE O
    O = [O,OUTPUT.O];
    
    
    % ------ LOCATE VF (Similar Procedure)
    ZZ = OUTPUT.Vf';
    ZZ = addVector(ZZ,year(i),SCALE,dt);
    Vf_val{pointi}(itag:itag+numel(ZZ)-1) = ZZ;
    
    % ------ LOCATE gwDisPart (Similar Procedure)
    ZZ = OUTPUT.gwDisPart;
    ZZ = addVector(ZZ,year(i),SCALE,dt);
    gwDisPart_val{pointi}(itag:itag+numel(ZZ)-1) = ZZ;
    
    % ------ LOCATE leakage (Similar Procedure)
    ZZ = OUTPUT.leakage;
    ZZ = addVector(ZZ,year(i),SCALE,dt);
    leakage_val{pointi}(itag:itag+numel(ZZ)-1) = ZZ;
    
    % ------ LOCATE ETr (Similar Procedure)
    ZZ = OUTPUT.ETr';
    ZZ = addVector(ZZ,year(i),SCALE,dt);
    ETr_val{pointi}(itag:itag+numel(ZZ)-1) = ZZ;
    
    
    itag = itag + numel(ZZ);
    
    fprintf('Load %2d successfully!\n',i);
    
end

% SAVE
save(goalFileName,'year','Runoff_time','Runoff_val','Qpoint_val','QpointS_val',...
    'Vf_val','gwDisPart_val','ETr_val','leakage_val','O');

end



%% AUXILLARY
function ZZ = addVector(ZZ,year,SCALE,dt)
% ADD MISSING
addZZ = NaN(sum(eomday(year,1:12))*24-length(ZZ),1);
ZZ = [ZZ;addZZ];

% (AGGREGATE)
if SCALE ~= dt
    ZZ = transpose(mean(buffer(ZZ(:),SCALE)));
end
end






