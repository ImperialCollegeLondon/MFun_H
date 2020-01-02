% CAT-PLOT-ARI-WITHR

function Cat_PLOT_ARI_withRObs_4Season(dt,CatchName,sp,ha,shift,SEASON,YRANGE)
% Example:
% sp = 0;%saveplot or not
% Cat_PLOT_ARI_withRObs_4Season(1,'Welland',sp)
% Update: 2019.12.29 Yuting

% CatchName = 'Welland';%'Stratford';%

fv = @(T)T;%-log(-log(1-1./T));%
a = 0.44;
getGringP = @(x)((1:length(x))-a)./(length(x)+1-2*a);
getReturnY = @(x) 1./(1-squeeze(x));

mermethod = 'GEAR';%'GEAR';%'BK';%'KED';%'CKED';%
[sim,obs] = Cat_PLOT_ARI(dt,{mermethod},CatchName);% only BK is considered for East Anglia

% prepare for Running R.scirpt
createSeasonTable=@(AM)table(AM(1,:)',AM(2,:)',AM(3,:)',AM(4,:)','VariableNames',{'MAM','JJA','SON','DJF'});
A = createSeasonTable(obs.AM_o);
writetable(A,['C:\Users\Yuting Chen\Desktop\Extreme\matlab_batch\obs_in_',...
    num2str(dt),'_',mermethod,'_',CatchName,'.txt']);

sim.AM_s = sim.AM_s;
sim_aux = permute(sim.AM_s,[1 3 2]);

sim_aa = [];
for season_i = 1:4
    sim_aa(:,season_i) = reshape(squeeze(sim_aux(:,:,season_i)),[],1);
end
sim_aux = sim_aa;

A = createSeasonTable(sim_aux');
writetable(A,['C:\Users\Yuting Chen\Desktop\Extreme\matlab_batch\sim_in_',...
    num2str(dt),'_',CatchName,'.txt']);
% close all;


%%
rfpath = 'C:\Users\Yuting Chen\Desktop\Extreme\matlab_batch\';
obsCol = [1 0 0];
simCol = [0.4 0.4 0.4];

for pltag = 1:length(SEASON)
    axes(ha(get_pos(pltag,shift)))

    %     figure;
    %     setFigureProperty;
    %     set(0,'defaultAxesFontSize',10);
    season = SEASON(pltag);
    %% PLOT SIM;
    % simulated series;
    sim_aux_x = [];
    sim_aux_y = [];
    
    
    for iter = 1:size(sim.AM_s,1)
        
        sim_aux_x(iter,:) = getReturnY(getGringP(sim.AM_s(iter,season,:)));
        sim_aux_y(iter,:) = sort(squeeze(sim.AM_s(iter,season,:)));
        
        [parmhat,~,~] = fitGEV(sim_aux_y(iter,:),'method','moments');
        rtx_temp = 1.1:0.1:100;
        y_temp = gevinv(1-1./rtx_temp,parmhat(1),parmhat(2),parmhat(3));
        
        % PLOT %
        if 0
        hsimfit = plot(fv(rtx_temp),y_temp,'-',...
            'Markersize',3,'color',simCol,'markerfacecolor',simCol,'linewidth',1);
        hold on;
        set(hsimfit,'visible','off')
        hsim = plot(fv(sim_aux_x(iter,:)),sim_aux_y(iter,:),'s',...
            'Markersize',2,'color',simCol,'markerfacecolor',simCol,'linewidth',1);
        hold on;
        set(hsim,'visible','off')
        end
    end
%     hsim = fillArea(fv(sim_aux_x(1,:)),nanmin(sim_aux_y),...
%         nanmax(sim_aux_y),simCol);
    % PLOT SIM line%
    sim_aux_x = getReturnY(getGringP(reshape(squeeze(sim.AM_s(:,season,:)),1,[])));
    sim_aux_y = sort(reshape(squeeze(sim.AM_s(:,season,:)),1,[]));

    hsim = plot(fv(sim_aux_x),sim_aux_y,'s',...
        'Markersize',3,'color',simCol,'markerfacecolor',simCol,'linewidth',1);
    set(hsim,'visible','off')
    hold on;
    
    xx = [1.1:0.1:100]';
 
    areacolor = simCol;
    rt = sim_aux_x;

    % load OBS %
    [rt,lmu,minS,maxS,areacolor] = loadOBS();
    
    % Plot OBS %
    [hobs,hobsGumbel,hobsGumbeluncer] = plotOBS(0.1);
    
    [rt,lmu,minS,maxS,areacolor] = loadSIM();
    hsimGumbel = plotSIM(0.1);
    
    % Legend %
    legend([hsimGumbel,hsimGumbeluncer,hobs,hobsGumbel,hobsGumbeluncer],...
        'Input: Sim Rain','95%CI','Input: Obs Rain','Gumbel Fit','95%CI',...
        'Location','Northwest');% hobs,'Observed',
    legend('boxoff')
    xlabel('Return Period [year]','FontSize',12);
    ylabel([timescale(dt),'Peak Discharge[m^3/s]'],'FontSize',12);
    ax = gca;
    hold on;
    grid minor
    % ax.XAxis.Scale = 'log';
    ylim([0,YRANGE(pltag)]);
    xlim([fv(1.1) fv(22)])
    
    text(fv(20),max(ylim)*0.1,[CatchName,'-',getSeasonName(season)],...
        'background',[0.9 0.9 0.9],...
        'HorizontalAlignment','right',...
        'VerticalAlignment','middle');
    
    XYWH = [150,150,400,200];
    set(gcf,'units','points','position',XYWH);
    
    
    % ax.XAxis.Scale = 'log';
end
set(ha(2:end),'YLabel',[]);

if sp
    reply = input('Do you want to save Y/N [Y]:','s');
    if strcmp(reply,'Y')
        
        savePath = ['C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\STNSRP\FIG_Pub\'];
        saveName = [savePath,'FFC_',CatchName,'_',num2str(dt),'h_24sim_',num2str(length(ha)),'Seasons_DGEAR'];
        saveas(gcf,saveName,'fig');
        saveas(gcf,saveName,'png');
        saveas(gcf,saveName,'eps');
        fprintf('Saved.');
        
    end
end



    function str = timescale(dt)
        switch(dt)
            case 1
                str = 'Hourly';
            case 24
                str = 'Daily';
        
        end
    end
    function hsimGumbel = plotSIM(alpha)
        hsimGumbel = plot(fv(xx),getGumbelYY_haveRT(rt,lmu(:,2),xx),'-',...
            'color',simCol,'linewidth',2);hold on;
        if 1
            hsimGumbeluncer = fillArea(fv(xx),minS,maxS,areacolor,alpha);
        elseif 1
            % not plot simulated CI because it doesn't realy be informative
            % for my aim..
            hsimGumbeluncer = plot(fv(xx),minS,'--','color',areacolor);hold on;
            plot(fv(xx),maxS,'--','color',areacolor);hold on;
        end
        text10_15(fv(10),fv(15),fv(1.1),lmu,simCol);
        
    end

    function [ul,ml,ll] = getprc(sim_aux_y)
        ul = prctile(sim_aux_y,95);
        ml = prctile(sim_aux_y,50);
        ll = prctile(sim_aux_y,5);
    end

    function [hobs,hobsGumbel,hobsGumbeluncer] = plotOBS(alpha)

        hobsGumbel = plot(fv(xx),getGumbelYY_haveRT(rt,lmu(:,2),xx),'-',...
            'color',obsCol,'linewidth',1);hold on;
        if 1
            hobsGumbeluncer = fillArea(fv(xx),minS,maxS,areacolor,alpha);
        else
            hobsGumbeluncer = plot(fv(xx),minS,'-.','color',areacolor);hold on;
            plot(fv(xx),maxS,'-.','color',areacolor);hold on;
        end
        % true p from observed
        hobs = plot(fv(getReturnY(getGringP(obs.AM_o(season,:)))),...
            sort(obs.AM_o(season,:)),'x',...
            'color',obsCol,'Markersize',3,'markerfacecolor',obsCol,'linewidth',1);
        text10_15(fv(10),fv(15),fv(1.1),lmu,obsCol);
        % set(hobs,'visible','off')
        % plot(grv(:,[1 3]),lmu(:,[1 3]),'b:','linewidth',1);hold on;
    end

    function [rt,lmu,minS,maxS,areacolor] = loadOBS();
        rt = [1.1:100.1]';
        rfName = ['obs_out_',num2str(dt),'_',CatchName,'_',mermethod,'_',getSeasonName(season),'.txt'];
        lmu = load([rfpath,rfName]);
        grv = repmat(rt,1,3);
        minS = getGumbelYY_haveRT(rt,lmu(:,1),xx);
        maxS = getGumbelYY_haveRT(rt,lmu(:,3),xx);
        areacolor = obsCol;
    end

    function [rt,lmu,minS,maxS,areacolor] = loadSIM();
        rt = [1.1:100.1]';
        rfName = ['sim_out_',num2str(dt),'_',CatchName,'_',getSeasonName(season),'.txt'];
        lmu = load([rfpath,rfName]);
        grv = repmat(rt,1,3);
        minS = getGumbelYY_haveRT(rt,lmu(:,1),xx);
        maxS = getGumbelYY_haveRT(rt,lmu(:,3),xx);
        areacolor = simCol;
    end
    function id = get_pos(season,shift)
    id = pltag*2-(2-shift)+1;
    end

end

function [parmhat,hobs,hfit] = getFitAMD(AM0,color,masize)

if isvector(AM0)
    
    [parmhat,hobs,hfit] = fitGEV(AM0,'dataPlot',1,'colorI',color,...
        'masize',masize,...
        'method','Gringorten');
else
    for i = 1:size(AM0,1)
        [parmhat,hobs,hfit] = fitGEV(AM0(i,:),'dataPlot',1,'colorI',color,...
            'method','Gringorten');
    end
end

end

function text10_15(x1,x2,inix,lmu,col)

% % add line;
% plot([x1,x1,inix],[0,lmu(9,2),lmu(9,2)],'k:','color',col);hold on;
% plot([x2,x2,inix],[0,lmu(14,2),lmu(14,2)],'k:','color',col);hold on;
% plot(x1,lmu(9,2),'ko','markerfacecolor',[0.5 0.5 0.5],'markersize',4);hold on;
% plot(x2,lmu(14,2),'ko','markerfacecolor',[0.5 0.5 0.5],'markersize',4);hold on;
% 
% % add text;
% rt = 10;
% text(x1,lmu(rt-1,2),sprintf('%.2f',lmu(rt-1,2)),'fontsize',12,...
%     'HorizontalAlignment','left',...
%     'VerticalAlignment','bottom',...
%     'Color',col);
% rt = 15;
% text(x2,lmu(rt-1,2),sprintf('%.2f',lmu(rt-1,2)),'fontsize',12,...
%     'HorizontalAlignment','left',...
%     'VerticalAlignment','bottom',...
%     'Color',col);

end

function h_fill = fillArea(x,minS,maxS,color,alpha)
%{
h_fill = plot(x,minS,'--','color',color);
plot(x,maxS,'--','color',color);
%}
if ~ischar(color)
    areaColor = color.*ones(1,3);
else
    areaColor = color;
end
p = plot(x,minS,x,maxS);
% ylim([10 60]);
YLIM = get(gca,'YLim');   delete(p);
h_fill = fill([x',x(end:-1:1)'],[maxS',minS(end:-1:1)'],areaColor,...
    'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',alpha);
hold on;

end

function yy = getSimYY(rt,y,rtxx)

op=1-1./rt;
x = -log(-log(op(:)));

opp = 1-1./rtxx;
xx = -log(-log(opp(:)));
yy = interp1(x,y,xx);

end

function yy = getGumbelYY_haveRT(rt,y0,rtxx)
op=1-1./rt;
p = polyfit(-log(-log(op(:))),sort(y0(:)),1);
opp = 1-1./rtxx;
yy = polyval(p,-log(-log(opp(:))));

end


function [sim,obs] = Cat_PLOT_ARI(dt,method,CatchName)
% dt = 24;
mon = 1:12;
annualPL = false;

% Go back; use original data

% method = {'BK'};
xx = [0.01:0.2:100];
[timeSeries] = Cat_PLOT_CompareWithNS_maxArea(method,dt,mon,annualPL,-1,xx,CatchName);

sim = struct;
obs = struct;

iter = length(timeSeries.sim_val);

% sim
for ii = 1:iter
    
    % Estimate empirical CDF
    C = timeSeries.sim_val{ii};
    am = findFlowSeasonMax(timeSeries.sim_time,C);
    for season = 1:size(am,2)
        [sim.AM_s(ii,season,:),sim.F_s(ii,season,:)] = ...
            getAM_F(reshape(am(:,season),1,[]));
    end
    
end

% radar
am = findFlowSeasonMax(timeSeries.radar_time, timeSeries.radar_val);
for season = 1:size(am,2)
    [obs.AM_o(season,:),obs.F_o(season,:)] = getAM_F(reshape(am(:,season),...
        1,[]));
end

end

function name = getSeasonName(season)
switch (season)
    case 1
        name = 'MAM';
    case 2
        name = 'JJA';
    case 3
        name = 'SON';
    case 4
        name = 'DJF';
end
end

function out = findFlowSeasonMax(ts,C)
ts = reshape(ts,[],1);
C = reshape(C,[],1);
try
    out = [];
    dateV = datevec(ts);
    yearRange = reshape(unique(dateV(:,1)),1,[]);
    for year = yearRange
        C1 = C(dateV(:,1) == year);
        T1 = ts(dateV(:,1) == year);
        out = [out;getSeasonMax_aux(T1,C1)];
    end
catch
    1;
end

    function out_aux = getSeasonMax_aux(T1,C1)
        d_aux = datevec(T1);
        mon = d_aux(:,2);
        season_ind = findSeasonalInd(mon);
        out_aux = cell2mat(cellfun(@(x)max(C1(x)),season_ind,...
            'UniformOutput', false));
    end

    function [season_ind] = findSeasonalInd(month)
        season_ind = cell(1,4);
        season_ind{1} = find(month==3 | month==4 | month==5);% MAM
        season_ind{2} = find(month==6 | month==7 | month==8);% JJA
        season_ind{3} = find(month==9 | month==10 | month==11);% SON
        season_ind{4} = find(month==12 | month==1 | month==2);% DJF
%         season_ind = cell(1,2);
%         season_ind{1} = find(month==6 | month==7 | month==8 | month == 9);% JJAS
%         season_ind{2} = find(month==12 | month==1 | month==2 | month == 3);% DJFM
    end
end


function [am_,f_] = getAM_F(am)
a = 0.44;
am_ = sort(am);
f_ = ((1:length(am))-a)/(length(am)+1-2*a);
end
