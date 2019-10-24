function [h] = Plot_RR(RUNOFF,RAINFALL,Unit)
% PLOT_RR plot the rainfall/runoff time series in one figure;
% Upper axes for rainfall
% Lower axes for runoff

% Input: RUNOFF: <cell> storing for runoff time series;
%        RAINFALL: <cell> storing rainfall time series;
%        Unit: <string>
%                 'Daily'
%                 'Hourly'
%
% Output:handle
%
% Example:
%        model = arima('Constant',0.01,'AR',{0.7,0.25},'Variance',.1);
%        rng('default')
%        N = 200;
%        Rain = simulate(model,N);
%        Runoff1 = Rain/10 + rand(N,1)/10;
%        Runoff2 = Rain([3:end,1:2])/10;
%        Plot_RR({Runoff1,Runoff2},Rain,'Daily');

% by Yuting Chen
%    yuting.chen17@imperial.ac.uk
%    Civil and Environmental Engineering
%    Imperial College London

if ~iscell(RUNOFF)
    RUNOFF = {RUNOFF};
end

if ~iscell(RAINFALL)
    RAINFALL = {RAINFALL};
end

setFigureProperty('Full');

RONum = length(RUNOFF);
SP = {'ko:','r-','b.:'};
if RONum > 0 && RONum < 3
    YM = 0;
    for i = 1:RONum
        
        rf = RUNOFF{i};
        
        h(i) = plot(rf,SP{i},'Linewidth',1,'Markersize',3);
        
        hold on;
        YM = max(YM,max(rf)*2.2);
    end
    
    xlim([1 length(rf)]);
    ylim([0 YM]);
    setYlabel(Unit,'Runoff');
elseif RONum >= 3
    YM = 0;
    for i = 1:RONum-1
        
        rf = RUNOFF{i};
        h(i) = plot(rf,'-','Linewidth',1,'Markersize',3);
        hold on;
        YM = max(YM,max(rf)*2.2);
        
    end
    rf = RUNOFF{RONum};
    h(RONum) = plot(rf,'ko:','Linewidth',1,'Markersize',2);
    hold on;
    YM = max(YM,max(rf)*2.2);
    xlim([1 length(rf)]);
    ylim([0 YM]);
    setYlabel(Unit,'Runoff');
else
    
    fprintf(['Check Input RUNOFF\n']);
    
    return;
    
end

color2 = [0 0.4470 0.7410];

ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
set(ax1,'box','off','color','none')
axes('Position',get(ax1,'Position'),'box','on','xtick',[],'ytick',[]);
axes(ax1);
ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','right','Color','none');
ax2.XColor = color2;
ax2.YColor = color2;
hold on

YM = 0;
if length(RAINFALL) == 1
    rf = RAINFALL{1};
    h(1+RONum) = bar(rf,'FaceColor',color2,'Parent',ax2);
    
    hold on;
    YM = max(YM,max(rf)*2.2);
else
    for i = 1:length(RAINFALL)
        
        rf = RAINFALL{i};
        h(i+RONum) = bar(rf,'Parent',ax2,'FaceAlpha', 0.7);
        
        hold on;
        YM = max(YM,max(rf)*2.2);
    end
end

ylim([0 YM])
xlim([1 max(length(rf),length(rf))]);

set(ax2,'Ydir','reverse')
setYlabel(Unit,'Rainfall');

grid minor

end



function setYlabel(unit,type)
switch(unit)
    case 'Daily'
        switch (type)
            case 'Rainfall'
                
                ylabel('Daily Rainfall[m^3/s]');
            case 'Runoff'
                ylabel('Gauge Daily Flow[m^3/s]');
                xlabel('Time [d]');
        end
    case 'Hourly'
        switch (type)
            case 'Rainfall'
                ylabel('Hourly Rainfall[mm/h]');
            case 'Runoff'
                ylabel('Gauge Hourly Flow[mm/h]');
                xlabel('Time [hour]');
        end
        
    otherwise
        
end

end

