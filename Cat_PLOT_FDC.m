CatchName =  'Welland';%'Stratford';%
close all

expNos = [201:220];
annualPL = false;
figure;
setFigureProperty('Paper')
ha = tight_subplot(1,2,[0.05 0.05],[.12 0.01],[.15 .10]);
pl = 1;
for seasonName = {'JJAS','DJFM'}%,'All'
% seasonName = 'JJAS';% seasonName = 'ALL';
axes(ha(pl))
seasonName = seasonName{1};
mon = getMons(seasonName);

itag = 1;
for dt = 24%[1,24]
    disFlow = 0.1;%0.2;
    
    method = {'GEAR'};%{'BK'};%{'CKED'};
    xx = [0.001:0.001:100];
    [timeSeries] = Cat_PLOT_CompareWithNS_maxArea(method,dt,mon,annualPL,disFlow,xx,CatchName,expNos,'FDC',true);
    
    itag = itag+1;
    if dt == 1
        ylabel('Hourly Flow (m^3)/s')
    elseif dt == 24
        ylabel('Daily Flow (m^3)/s');
    end
    xlabel('Percentage of time flow exceeded');
    text(2,0.2,...
        [char(getMonthName(mon(1))),' to ',char(getMonthName(mon(end)))],...
        'fontsize',8);
    set(gca,'linewidth',1)
    xlim([0 100]);
    ylim([disFlow 600]);
    
    yticks([-0.1,1,10,100])
    yticklabels([-0.1,1,10,100])
    xticks(0:20:100)
    xticklabels([0,20,40,60,80,100])
    axis('square')
    box off
    pl = pl+1;
end

end
format_xylabel(ha,1,2)
filePath = ['C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\STNSRP\FIG_Pub'];
% filename = [filePath,filesep,'FDC_',CatchName,'Area_',num2str(dt),'h_20sim_',...
%     seasonName,'_',method{1}];
filename = [filePath,filesep,'FDC_',CatchName,'Area_',num2str(dt),'h_20sim_',method{1}];
savePlot(filename,'targetSize','1c','needreply','Y');


function mon = getMons(seasonName)
switch(upper(seasonName))
    case 'DJFM'
        mon = [12 1 2 3];
    case 'JJAS'
        mon = 6:9;
    case 'ALL'
        mon = 1:12;
end
end

