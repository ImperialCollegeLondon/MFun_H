CatchName = 'Welland';%'Stratford';%

mon = [12 1 2 3];
seasonName = 'DJFM';

% mon = [6:9];
% seasonName = 'JJAS';

annualPL = false;
itag = 1;
for dt = 1%24%[1,24]%[24]
    figure;
    setFigureProperty;
    disFlow = 0.2;
    
    method = {'GEAR'};%'GEAR'{'BK'};%{'CKED'};
    xx = [0.001:0.001:100];%[0.001:0.001:100];
    [timeSeries] = Cat_PLOT_CompareWithNS_maxArea(method,dt,mon,annualPL,disFlow,xx,CatchName);
    
    itag = itag+1;
    if dt == 1
        ylabel('Hourly Flow (m^3)/s')
    elseif dt == 24
        ylabel('Daily Flow (m^3)/s');
    end
    xlabel('Percentage of time flow exceeded');
    text(2,0.3,...
        [char(getMonthName(mon(1))),' to ',char(getMonthName(mon(end)))],...
        'fontsize',10);
    
    xlim([0 100]);
    ylim([disFlow 600]);
    XYWH = [150,150,250,225];
    set(gcf,'units','points','position',XYWH);
    savetag = input('Ensure to save this figure? Y/N','s');
    
    if strcmp(savetag,'Y')
        
        savePath = ['C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\STNSRP\FIG_Pub\'];
        saveName = ['FDC_',CatchName,'Area_',num2str(dt),'h_10sim_',...
            seasonName,'_',method{1}];
        saveas(gcf,[savePath,saveName],'fig');
        saveas(gcf,[savePath,saveName],'png');
        saveas(gcf,[savePath,saveName],'eps');
        fprintf('Saved');
        
    else
        fprintf('No Save');
    end
    
end


