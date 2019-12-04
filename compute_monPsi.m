function [mon_psi,xx,yy] = compute_monPsi(A,CompleteTS,GEAR,TOPOFILE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All info needed are saved in GEAR.mat, according to the different
% calendar month.
for MON = 1:12
    
    % contourf(GEAR(MON).XX,GEAR(MON).YY,GEAR(MON).RAIN);
    dt = 1;
    param = A.Result0_spatial(MON,:);
    psi_val = param(8:end-2);
    psi_x = CompleteTS.ENA(:,1)';
    psi_y = CompleteTS.ENA(:,2)';
    
    % FIND corresponding fitted mean rainfall using NS(s-t) in continuous version
    %     IF1 = @(z)z.*((1-z).*param(2)+z.*param(3))./param(5);
    %     rain_val = param(1)*psi_val*1*quadgk(IF1,0,1)*24*365;
    
    [rain_val, ~, ~, ~, ~, ~, ~, ~, ~] = ns_continuousStormType_model_statistics_exp( param(1), param(2), ...
        param(3), param(4),param(5), param(6), param(7), psi_val, dt, [], 'NoPDD & NoCcorr');
    rain_val = rain_val * 24 * 365; % 365 days are used due to gear_val is for 365 days
    
    gearX_aux = GEAR(MON).XX/1000;
    gearY_aux = GEAR(MON).YY/1000;
    gear_val = GEAR(MON).RAIN;
    pl = 0;
    [mat_psi] = getScalingMatrix(pl,gearX_aux,gearY_aux,gear_val,psi_x,psi_y,psi_val,rain_val);
    
    %%
    % DETERMINE DOMAIN RANGE for CATCHMENT RESEARCH
    DATA = load(TOPOFILE,'x_yr','y_yr');
    CAT_x =  unique(round(DATA.x_yr/1000));
    CAT_y = unique(round(DATA.y_yr/1000));
    [xx,yy] = meshgrid(CAT_x,CAT_y(end:-1:1));
    CAT_psi = 0*CAT_x;
    for i = 1:size(xx,1)
        for j = 1:size(xx,2)
            dist_aux = sqrt((gearX_aux-xx(i,j)).^2+(gearY_aux-yy(i,j)).^2);
            ii_aux = find(dist_aux == min(dist_aux(:)));
            CAT_psi(i,j) = mat_psi(ii_aux(1));
        end
    end
    
    mon_psi(MON,:) = CAT_psi(:)';% CAT_psi(~isnan(DATA.DEM_yr));
    % contourf(CAT_psi.*mean(rain_val./psi_val))
    
end



end