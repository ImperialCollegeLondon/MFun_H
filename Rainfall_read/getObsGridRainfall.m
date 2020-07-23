function Rain_val = getObsGridRainfall(loc,time,dataSource)

xx = loc.x;
yy = loc.y;

if isfield(time,'month')
    month = time.month;
end

switch(dataSource)
    
    case 'GEAR'
        GEAR = getfield(load('GEAR.mat'),'GEAR');% mean monthly in [2007,2015];
        xyR = GEAR(month);
        
        xy_ul = [xyR.XX(1,1),xyR.YY(1,1)];
        dx = xyR.XX(1,2)-xyR.XX(1,1);
        dy = xyR.YY(2,1)-xyR.YY(1,1);
        
        ii = ceil((yy-xy_ul(2))/dy+0.5);
        jj = ceil((xx-xy_ul(1))/dx+0.5);
        try
            Rain_val = xyR.RAIN(ii(:,1),jj(1,:));
        catch
            1;
        end
    otherwise
        % ...
        % ... maybe NIMROD radar or other data source here.
        % ...
end

end