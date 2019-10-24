function PET = Interpolate2D_timeseries(PET,MI)
%% INTERPOLATE2D_TIMESERIES finish the interpolation for time series over the matrix.

% Missing data is shown as NaN in the matrix PET
% 
% Input: PET: 3D Matrix<double> 
%             X * Y * T
%        MI:  Vector<double>
%
% Output:PET: 3D Matrix<double> without NaN value;
% 
% Method: Nearest neighbor interpolation in GRIDDATAN
%
% by Yuting Chen
%    Imperial College London
%    yuting.cehn17@imperial.ac.uk

fprintf('Interpolate2D_timeseries started....\n');

if isempty(MI)
    fprintf('Check: MI is empty\n (still maybe right)');
    return
end
if MI(1) == 1
    error('Please check PET data, in this case, interpolation is not processed.');
end

% PET_ori = PET;

iForFill = [];
for i = 1:length(MI)
    disp(i);
    if i ~= length(MI)
        iForFill = [iForFill,MI(i)];
        if MI(i+1) == MI(i)+1
            %%%
        else
            % Interpolate between iForFill(1)-1, iForFill(end)+1
            RI = iForFill(1)-1:iForFill(end)+1;
            [xx,yy,zz] = meshgrid(1:size(PET,2),1:size(PET,1),RI);
            X = [xx(:),yy(:),zz(:)];
            v = PET(:,:,RI);
            v = v(~isnan(v));
            X = X(~isnan(v),:);
            [xx,yy,zz] = meshgrid(1:size(PET,2),1:size(PET,1),iForFill);
            xq = [xx(:) yy(:) zz(:)];
            vq = reshape(griddatan(X,v,xq,'nearest'),size(xx));
            PET(:,:,iForFill) = vq;
            
            %             vq = reshape(vq,size(xx));
            %             plot3(X(:,1),X(:,2),X(:,3),'.')
            %             hold on
            %             slice(xx,yy,zz,vq,[0.2 0.4 0.6 0.8],0.5,0.5)
            
            iForFill = [];
            
            if length(RI)>10
                error('1');
            end
        end
    else
        iForFill = [iForFill,MI(i)];
        for rii = 1:length(iForFill)
            PET(:,:,iForFill(rii)) = PET(:,:,iForFill(1)-1);
        end
        
    end
    
    
end

fprintf('Interpolate2D_timeseries finished.\n');

end