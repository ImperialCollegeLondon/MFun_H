CatchName ='Stratford';%'Welland';% 

mon = [12 1 2 3];
seasonName = 'DJFM';

% mon = [6:9];
% seasonName = 'JJAS';

% mon = [1:12];
% seasonName = 'ALL';
close all

annualPL = false;
itag = 1;
for dt = 1%[1,24]
    figure;
    setFigureProperty('Subplot2')
    disFlow = 0.1;%0.2;
    
    method = {'GEAR'};%{'BK'};%{'CKED'};
    xx = [0.001:0.001:100];%[0.001:0.001:100];
    [timeSeries] = Cat_PLOT_CompareWithNS_maxArea(method,dt,mon,annualPL,disFlow,xx,CatchName);
    
    itag = itag+1;
    if dt == 1
        ylabel('Hourly Flow (m^3)/s')
    elseif dt == 24
        ylabel('Daily Flow (m^3)/s');
    end
    xlabel('Percentage of time flow exceeded');
    text(2,0.2,...
        [char(getMonthName(mon(1))),' to ',char(getMonthName(mon(end)))],...
        'fontsize',12);
    set(gca,'linewidth',2)
    xlim([0 100]);
    ylim([disFlow 600]);
    
    yticks([-0.1,1,10,100])
    yticklabels([-0.1,1,10,100])
    xticks(0:20:100)
    xticklabels([0,20,40,60,80,100])
    
    filePath = ['C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\STNSRP\FIG_Pub'];
    filename = [filePath,filesep,'FDC_',CatchName,'Area_',num2str(dt),'h_10sim_',...
        seasonName,'_',method{1}];
    savePlot(filename,'XYWH',[150,150,250,225],'needreply','N');
    
end


