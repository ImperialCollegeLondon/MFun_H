function [mat_psi] = getScalingMatrix(pl,gearX,gearY,gear_val,psi_x,psi_y,psi_val,rain_val)
% GETSCALINGMATRIX return a matrix of scaling factor for NS(s-t) based model
% Reference data: GEAR monthly radar estimate
%                 <....link....>
% Input: pl
%        gearX
%        gearY
%        gear_val: mm/365 days
%        psi_x
%        psi_y
%        psi_val
%        rain_val
% OUTPUT:mat_psi
%
% Example:.....
%         .....
%
% Description: IDW method are used to determine the scaling factor for each
% grid
% Yuting CHEN
% Update: 11/03/2019
% Note: Here, the point statistics at each stations are remained.


% FIND corresponding index of PSI- vector in GEAR- matrix
ind = NaN(1,numel(psi_x));
for i = 1:numel(psi_x)
    x_aux = psi_x(i);
    y_aux = psi_y(i);
    dis_2_aux = (gearX(:)-x_aux).^2+(gearY(:)-y_aux).^2;
    ind_aux = find(dis_2_aux == min(dis_2_aux));
    ind(i) = ind_aux(1);
end

% GENERATE PSI Matrix
% method: IDW
mat_psi = 0*gear_val;
for i = 1:size(gearX,1)
    for j = 1:size(gearX,2)
        w_aux = 1./sqrt((psi_x-gearX(i,j)).^2+(psi_y-gearY(i,j)).^2);
        w_aux = w_aux./sum(w_aux);
        mat_psi(i,j) = sum((gear_val(i,j)./(rain_val./psi_val)).*w_aux);
    end
end

mm = mean(rain_val./psi_val);

if pl
    figure;
    setFigureProperty;
    subplot(131)
    contourf(mat_psi*mm);colorbar
    subplot(132)
    contourf(gear_val);colorbar
    subplot(133);
    contourf(mat_psi*mm-gear_val);colorbar
end

end
