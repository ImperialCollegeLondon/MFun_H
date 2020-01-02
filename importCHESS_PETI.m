%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% download CHESS data & Import some CHESS PETI;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% https://catalogue.ceh.ac.uk/datastore/eidchub/8baf805d-39ce-4dac-b224-c926ada353b7/

function [PETGroup,status] = importCHESS_PETI(XX,YY,YEAR)
% Example:
%     %Import PETI data at Yare station
%     XX = 618200;
%     YY = 308200;
%     YEAR = 2007:2008;
%     [PET,status] = importCHESS_PETI(pointX,pointY,YEAR)
%     status

% Output:
%     PETGroup: <3d double vector>
%          [mm/day]
%          [daily PET i]
%     status: 1: successfully output;
%             0: fail

%------------------------------%
% by Yuting Chen                  %
% Imperial College London      %
% yuting.chen17@imperial.ac.uk %
%------------------------------%

%%%%% DOWNLOAD+SAVE CHESS DATA;
% clear;clc
% % mkdir([cd,'\CHESS_PETI']);
% %%
% YEARRANGE = [2007:2014];
% options = weboptions('username','yutingchen0604@hotmail.com','password','AaBb14207','Timeout',Inf);
% weblink = 'https://catalogue.ceh.ac.uk/datastore/eidchub/8baf805d-39ce-4dac-b224-c926ada353b7/peti/';
% for year = YEARRANGE
%     for mon = 1:12
%         try
%             OFN = websave(['C:\Users\Yuting Chen\Box\CHESS_PETI\chess_peti_wwg_',num2str(year),sprintf('%02d', mon),'.nc'],...
%                 [weblink,'chess_peti_wwg_',num2str(year),sprintf('%02d', mon),'.nc'],options);
%             disp(['Finished:',num2str(year),'-',sprintf('%02d', mon)]);
%         catch
%             disp('111');
%         end
%     end
% end

%% generate CHESS PETI [mm/day] data for *** region;
status = 1;

fprintf('importCHESS_PET started....');

PETGroup = NaN(size(XX,1),size(XX,2),2);

findSDi = @(X,px)find(abs(X-px) == min(abs(X-px)));

auxTag = 1;
try
    for year = YEAR
        for mon = 1:12
            source = ['K:\CHESS_PETI\chess_peti_wwg_',num2str(year),sprintf('%02d', mon),'.nc'];
            % finfo = ncinfo(source);
            peti = ncread(source,'peti');
            x = ncread(source,'x');
            y = ncread(source,'y');
            for pi = 1:size(XX,1)
                for pj = 1:size(XX,2)
                    ix = findSDi(x,XX(pi,pj));
                    iy = findSDi(y,YY(pi,pj));
                    petAux = squeeze(mean(mean(peti(ix,iy,:),1),2));
                    PETGroup(pi,pj,(auxTag:auxTag+eomday(year,mon)-1)) = petAux;
                end
            end
            auxTag = auxTag + eomday(year,mon);
        end
    end
catch
    status = 0;
end

fprintf('importCHESS_PET finished.\n');

end
