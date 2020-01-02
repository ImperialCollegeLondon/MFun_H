%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% download CHESS data & Import some CHESS PET;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% https://catalogue.ceh.ac.uk/datastore/eidchub/8baf805d-39ce-4dac-b224-c926ada353b7/

function [PETGroup,status] = importCHESS_PET(XX,YY,YEAR)
% Example:
%     %Import PET data at Yare station
%     XX = 618200;
%     YY = 308200;
%     YEAR = 2007:2008;
%     [PET,status] = importCHESS_PET(pointX,pointY,YEAR)
%     status

% Output:
%     PETGroup: <3D double matrix>
%          [mm/day]
%          [daily PET i]
%     status: 1: successfully output;
%             0: fail

%------------------------------%
% by Yuting Chen                  %
% Imperial College London      %
% yuting.chen17@imperial.ac.uk %
%------------------------------%

% %%%% DOWNLOAD+SAVE CHESS DATA;
% clear;clc
% % mkdir([cd,'\CHESS_PET']);
% %%
% YEARRANGE = [2007:2014];
% options = weboptions('username','yutingchen0604@hotmail.com','password','AaBb14207','Timeout',Inf);
% weblink = 'https://catalogue.ceh.ac.uk/datastore/eidchub/8baf805d-39ce-4dac-b224-c926ada353b7/pet/';
% for year = YEARRANGE
%     for mon = 1:12
%         try
%             OFN = websave(['C:\Users\Yuting Chen\Box\CHESS_PET\chess_pet_wwg_',num2str(year),sprintf('%02d', mon),'.nc'],...
%                 [weblink,'chess_pet_wwg_',num2str(year),sprintf('%02d', mon),'.nc'],options);
%             disp(['Finished:',num2str(year),'-',sprintf('%02d', mon)]);
%         catch
%             disp('111');
%         end
%     end
% end
% 
%% generate CHESS PET [mm/day] data for *** region;
status = 1;

fprintf('importCHESS_PET started....');

PETGroup = NaN(size(XX,1),size(XX,2),2);

findSDi = @(X,px)find(abs(X-px) == min(abs(X-px)));

auxTag = 1;
try
    for year = YEAR
        for mon = 1:12
            %C:\Users\Yuting Chen\Box
            source = ['K:\CHESS_PET\chess_pet_wwg_',num2str(year),sprintf('%02d', mon),'.nc'];
            % finfo = ncinfo(source);
            pet = ncread(source,'pet');
            x = ncread(source,'x');
            y = ncread(source,'y');
            for pi = 1:size(XX,1)
                for pj = 1:size(XX,2)
                    ix = findSDi(x,XX(pi,pj));
                    iy = findSDi(y,YY(pi,pj));
                    petAux = squeeze(pet(ix(1),iy(1),:));
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



