
close all
clear all

load PLM_TOPO_DATA2.mat;

% SN = 0*SN;

load PLM_data_h;

DTM = DEM_pl2;
[R, S] = dem_flow(DTM,cellsize,cellsize);

% %%% Precipitation
Prt = zeros(5*24,1);
Prt(5:10) = 10;
Prt(48:53) = 10;

Ndt=length(Prt);
x=x_pl;
y=y_pl;
%%%%%%%%%%%%%%%%%%%%%%
[m,n] = size(DTM);
MASK=ones(size(DTM));
MASK(isnan(DTM))=0;
SLO = S;
clear S %% fraction

% SLO(SLO<0.0001)=0.0001;
SLO(SLO<0.01)=0.01;

SLO=SLO.*MASK;
% SN = zeros(size(DTM)).*MASK;
SN(isnan(SN))=0;
%SLO(SN==1)=0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Matrice Zs [mm]

Zs=zeros(m,n);
Zs(:,:) = 400.*MASK; %%  Soil Active depth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Osat =  (zeros(m,n) + 0.35).*MASK;
Ohy =   (zeros(m,n) + 0.1).*MASK;
Oel =   (zeros(m,n) + 0.11).*MASK;
%%%%%%%%%%%%%%%%%%
nvg= 2.8;
% nvg= 1.05;
avg = -6.3; %% [1/m]
Ks=  (zeros(m,n) + 20).*MASK; %%% [mm/h]
mvg = (zeros(m,n) + (1 - 1/nvg)).*MASK;
Kbot=  (zeros(m,n) + 0.01).*MASK; %%% [mm/h]
%%%%%%%%%%%%%%%%%%%%%%%%%%
kF = 10000;%   %%% [h]  acquifer constant
% mF= 100000; %320; %%% [mm]
mF= 100000; %320; %%% [mm]

f=(Osat-Ohy)./mF; %%% [1/mm]
%%%%%
OPT_UNSAT = 1; %%% Option unsaturated 1 -- saturated bottom 0
CF=1;%0; %% Acquifer Yes - No
aR=100; %%[-] anisotropy ratio
aT= 1000*(cellsize^2)/cellsize; %%[mm] Area/Contour length ratio
Area= (cellsize^2)*sum(sum(MASK)); %% [m^2]

%%%%% dth --- >>>>
dt_in_d = 3600;

dt = 1*60;%%%[s]
dth= dt/3600;% [h]
dti= 10; %%[s] Internal Time step for Surface Channel-flow Routing
dti2= 5; %%[s] Internal Time step for Surface Overland-flow Routing
%%%%%%%%%%%%%%%%%%%%%%%%
DEM = GRIDobj(x_pl(1,:),y_pl(:,1),DEM_pl);

FD  = FLOWobj(DEM);
FA = flowacc(FD);

WC = zeros(m,n); %%% [m]  width channel
WC(SN==1)=0.0015*sqrt(FA.Z(SN==1)*cellsize^2); %% [m]

WC=WC.*MASK;
NMAN_C=SN*0.03;
NMAN_H=0.10;
NMAN_C=NMAN_C.*MASK;%%[s/(m^1/3)] manning coefficient
NMAN_H=NMAN_H.*MASK;
MRO = 0.0; %0.002; %%[m]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vmax = zeros(m,n);
Vmax =Zs.*(Osat-Ohy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Qout=zeros(1,Ndt);
QoutS=zeros(1,Ndt);
Qpoint=zeros(1,Ndt);
QpointS=zeros(1,Ndt);
QoutC=zeros(1,Ndt);
QpointC=zeros(1,Ndt);
Upoint=zeros(1,Ndt);
UCh_point=zeros(1,Ndt);
S_ave=zeros(1,Ndt);
Oave=zeros(1,Ndt);
Ostd=zeros(1,Ndt);
Cov_OQ=zeros(1,Ndt);
Cov_OE=zeros(1,Ndt);
ENT=zeros(1,Ndt);
ENT2=zeros(1,Ndt);
FEN=zeros(1,Ndt);
EPE=zeros(1,Ndt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Outlet
Xout=205; Yout=277;

SLO(Yout,Xout)=0.02;

%%% Evapotranspiration
ETr = zeros(1,Ndt);
ETPt =  zeros(1,Ndt);
%%% Initial Condition
QsurR=zeros(m,n);
QchR=zeros(m,n);
%%%%%%%%%%%%
Vf=zeros(1,Ndt);
Vf(1)=0;
%%%%%%% Intial Oi
%Oi=(Osat-0.15).*MASK;   O = Oi;
Oi=0.15.*MASK;
O=Oi;
%%%%%%%%%%%%%%%%%%%
CK1=zeros(1,Ndt);
CK2=zeros(1,Ndt);

col='og';
%%%

tic ;
eff_soil_sat=zeros(size(DTM));

DTM_nml = numel(DTM);
[mm,nn]=size(DTM);
XX = -T + speye(mm*nn,nn*mm);

tH = 0*DTM;
t = 0*DTM;
t1 = 0;
t2 = 0;
SNl = logical(SN);
Y = zeros(size(DTM));

for i=2:Ndt
    
    disp(i)
    
    if i == 2
        
        dt_q = dt;
        dth= dt_q/3600;% [h]
        
        sub_n = dt_in_d/dt_q;
    else
        
        dt_q = 0.1*nanmin(3600*V(:)./Qi(:));
%         dt_q = 0.1*(nanmin(nanmin((1./(To./Zs./(cellsize*1000)))))*3600);     
        dt_q = nanmax(dt_q,dt);
        dt_q = dt_in_d/ceil(dt_in_d/dt_q);
%         
        
%         dt_q = dt;
        dth= dt_q/3600;% [h]
        
        sub_n = dt_in_d/dt_q;
%         aaa(i) = dth;
    end
    
    for sub_i = 1:sub_n
        
        
%         t1 = nanmin(nanmin(tH(SN==0)));
%         t2 = nanmin(nanmin(t(SN==1)));
        
        dti  = nanmax(t1/10,5);
        dti2 = nanmax(t2/10,5);
        
        dti = dt_q/ceil(dt_q/dti);
        dti2 = dt_q/ceil(dt_q/dti2);
        
        % Matrix Flow
        Pr= Prt(i)*MASK + QsurR; %%% [mm/h]  Precipitation
        ET= ETPt(i)*MASK; %%% [mm/h]  evapotranspiration
        %BETA=ones(m,n);
        BETA=(O(:,:)-Ohy)./(Oel-Ohy);  BETA(BETA>1)=1; BETA(BETA<0)=0;
        BETA(isnan(BETA))=0; %%%%% evapotranspiration factor
        ET=ET.*BETA; %[mm/h]  evapotranspiration
        
        ETr(i) =  sum(sum(ET))*(cellsize^2/Area); %[mm/h]  evapotranspiration
        Rh = Pr - Ks;  Rh(Rh<0)=0; %%% [mm/h] Horton Runoff
        
        If = Pr;
        If(If>Ks)=Ks(If>Ks);
        If=If*dth; %% [mm] Infiltration
        
        %%%%%%% First Layer %%%%
        VI= (O-Ohy).*Zs; %%% Volume First layer present [mm]
        V= VI + If;  %% Volume first update [mm]
        O = V./Zs + Ohy;  O(isnan(O))=0; %%
        O(O>Osat)=Osat(O>Osat); %% Soil Moisture first update []
        
        Ko = Ks.*(((O-Ohy)./(Osat-Ohy)).^(1/2)).*...
            (1-(1-((O-Ohy)./(Osat-Ohy)).^(1./mvg)).^mvg).^2; %%%%%% Unsaturated conductivity surface [mm/h]
        Ko(isnan(Ko))=0; %%
        if OPT_UNSAT == 1
            tS = (Osat).*1000*cellsize./Ko; %%% [h]
            tS(isnan(tS))=0; %%
        else
            tS = (Osat).*1000*cellsize./Ks; %%% [h]
            tS(isnan(tS))=0; %%
        end
        %kdtS= dth./tS; %%[]
        %kdtS(isnan(kdtS))=0; %%
        %kdtS(kdtS>1)=1;
        kdtS=1;
        
        if OPT_UNSAT == 1
            To = (aR.*Ko./f).*(1 -exp(-f.*Zs)); %%% Trasmsissivity [mm^2/h]
            To(isnan(To))=0;
        else
            Z= V./(Osat-Ohy);%%% [mm]
            Z(isnan(Z))=0;
            To= Ks.*(Z); %% [mm^2/h]
        end
        Ks_Z = Ks.*exp(-f.*Zs); Ks_Z(Ks_Z<Ks/10)= Ks(Ks_Z<Ks/10)/10; %% [mm/h] Saturated Conductivity Bottom
        Ks_b =  exp((log(Kbot) + log(Ks_Z))/2);%  [mm/h] Conductivity between layer n and bedrock
        if OPT_UNSAT == 1
            Lk= CF*Ks_b.*(((O-Ohy)./(Osat-Ohy)).^(1/2)).*...
                (1-(1-((O-Ohy)./(Osat-Ohy)).^(1./mvg)).^mvg).^2; %% Unsaturated conductivity surface [mm/h]
            Lk(isnan(Lk))=0;
        else
            Lk=CF*Ks_b; %%%  [mm/h]
            Lk(isnan(Lk))=0;
        end
        
        Qi = kdtS.*To.*SLO/aT; %%% [mm/h] %%% Lateral Outflow
        V= V-Qi*dth - Lk*dth- ET*dth; %%% Volume second update [mm]
        
        Rd = V-Vmax; Rd(Rd<0)=0; % Dunne Runoff [mm]
        V=V-Rd; EX = -V.*(V<0); V(V<0)=0; %%% Volume Third update
        Rd=Rd/dth; %%% Dunne Runoff [mm/h]
        
        %%%% Acquifer
        Vf(i) = Vf(i-1) + sum(sum(Lk*dth)) -sum(sum(EX)); %%% [mm] Total
        if Vf(i)<0
            disp('Time Step too Long')
            break
        end
        
        Qsur= Rd + Rh; %%% Surface flow [mm/h]
        
        CK1(i)= sum(sum(VI))/1000*(cellsize^2)  -  sum(sum(V))/1000*(cellsize^2) + sum(sum(Pr*dth))/1000*(cellsize^2) ...
            - sum(sum(ET*dth))/1000*(cellsize^2) + sum(sum(EX))/1000*(cellsize^2) ...
            -  sum(sum(Qsur*dth))/1000*(cellsize^2) -  sum(sum(Qi*dth))/1000*(cellsize^2) - sum(sum(Lk*dth))/1000*(cellsize^2);  %%[m^3]
        CK1(i)=1000*CK1(i)/Area; %%[mm]
        
        %%%% Subsurface Routing
        Qseep=Qi*dth.*(SN); %%[mm]
        Qi=Qi.*(1-SN); %% [mm/h]
        %%% SUBSURFACE VELOCITY
        %     [QiM]=Flow_Routing_Step2(DTM,T,Qi*dth); %%[mm]
        
        Qsurl = reshape(Qi*dth,DTM_nml,1);
        QiM =(XX)*Qsurl;
        QiM  = reshape(QiM , [mm nn]);
        
        QoutS(i)= (sum(sum(Qi*dth))- sum(sum(QiM))); %%%  %%[mm]
        QpointS(i)=  Qi(Yout,Xout)*dth;%%%  %%[mm]
        
        Qi_fall = 0*MASK;
        tH_store = 0*MASK;
        %%%%% OVERLAND FLOW
        Qsur = Qsur*dth; %%[mm]
        
        if sum(sum(Qsur))>0
            for jjj=1:1:dt_q/dti
                %%%% SURFACE VELOCITY
                %WC= WC0 + 2*Y*cot(alpha); %%% [m]
                %t=(cellsize.*(WC.^0.4).*(NMAN.^0.6))./(((cellsize^2)*Qsur/3600/1000).^0.4.*(SLO).^0.3); %%% [s]
                YH = (Qsur*dth/1000);% .*cellsize./WC ; %% [m]
                %%%%% Ponding - Microroughness
                YH = YH - MRO ; YH(YH<0)=0;
                %%%%
                tH= (cellsize.*NMAN_H)./(YH.^(2/3).*SLO.^0.5); %%% [s]
                tH(isnan(tH))=0;
                tH_store = tH_store + tH*dti/dt_q;
                kdt= dti./tH; %[]
                kdt(kdt>1)=1;
                %kdt(kdt>0)=1; %%% option no-routing
                kdt(isnan(kdt))=0;
                %%%% Surface Routing
                
                Qsurl = reshape(kdt.*Qsur,DTM_nml,1);
                QsurM =(XX)*Qsurl;
                
                QsurM  = reshape(QsurM , [mm nn]);
                
                QsurM2 = QsurM.*(1-SN); %%[mm]
                Qi_fall = Qi_fall + QsurM.*(SN); %% [mm]
                I_fall = sum(sum(QsurM.*(SN)));
                QsurM = QsurM2;
                QsurR = QsurM + (Qsur - kdt.*Qsur) ; %%[mm]
                
                Qout(i)= Qout(i) + (sum(sum(Qsur))- sum(sum(QsurR))-I_fall); %%%%% -- [mm]
                Qpoint(i)= Qpoint(i) + kdt(Yout,Xout)*Qsur(Yout,Xout); %%%%% -- [mm]
                
                Qsur=QsurR; %%[mm]
            end
        else
            YH =Qsur ; %% [m]
            tH= Inf*Qsur; %%% [s]
            Qout(i)= 0; %%%%% -- [mm]
            Qpoint(i)=0; %%%%% -- [mm]
            
            QsurR=0; %%[mm]
            Qi_fall=0; %%[mm]
        end
        
        Qch = QchR*dth + Qseep + Qi_fall; %%% [mm]
        
        t_store = 0*MASK;
        %%%
        
        if sum(sum(Qch))>0
            for jjj=1:1:dt_q/dti2
                %%%%%% SURFACE VELOCITY %%%%%%%%%%%
                Y(SNl) = (Qch(SNl)/1000).*cellsize./WC(SNl) ; %% [m]
                t(SNl)= (cellsize.*NMAN_C(SNl))./(Y(SNl).^(2/3).*SLO(SNl).^0.5); %%% [s]
                t(isnan(t))=0;
                t_store = t_store + t*dti2/dt_q;
                kdt= dti2./t; %[]
                kdt(kdt>1)=1;
                kdt(isnan(kdt))=0;
                
                Qchl = reshape(kdt.*Qch,DTM_nml,1);
                QchM=(XX)*Qchl;
                QchM = reshape(QchM, [mm nn]);
                QchR = QchM + (Qch - kdt.*Qch) ; %% [mm]
                
                QoutC(i)= QoutC(i) + (sum(sum(Qch))- sum(sum(QchR))); %%%%% -- [mm]
                QpointC(i)= QpointC(i) + kdt(Yout,Xout)*Qch(Yout,Xout); %%%%% -- [mm]
                Qch=QchR; %%[mm]
            end
        else
            Y =Qch ; %% [m]
            t= Qch*Inf; %%% [s]
            QoutC(i)= 0; %%%%% -- [mm]
            QpointC(i)=0; %%%%% -- [mm]
            QchR=0; %%[mm]
        end
        
        %% addaptive time step addition
        
        t1 = nanmin(nanmin(tH(SN==0)));
        t2 = nanmin(nanmin(t(SN==1)));
%         
%         dti  = nanmax(t1/10,30);
%         dti2 = nanmax(t2/10,2);
%         
%         dti = dt_q/ceil(dt_q/dti);
%         dti2 = dt_q/ceil(dt_q/dti2);
        %%
        
        QsurR = QsurR/dth; %%% [mm/h]
        QchR = QchR/dth; %%% [mm/h]
        QiM = QiM/dth; %%% [mm/h]
        
        V = V + QiM*dth  ; %%% Volume fourth update [mm]
        Rd2 = V-Vmax; Rd2(Rd2<0)=0; % Dunne Runoff 2 [mm]
        V=V-Rd2; %%% Volume fifth update
        O= V./Zs + Ohy; O(isnan(O))=0; %% % Soil Moisture update []
        QsurR= QsurR + Rd2/dth;  %%% [mm/h]
        
        %%%% Acquifer Routing
        Qsub= Vf(i)/kF; %%% [mm/h]
        Vf(i)=Vf(i)- Qsub*dth; %%[mm]
        QchR= QchR + SN*Qsub/sum(sum(SN));  %%% [mm/h]
        
        CK2(i)= sum(sum(VI))/1000*(cellsize^2)  -  sum(sum(V))/1000*(cellsize^2) + sum(sum(Pr*dth))/1000*(cellsize^2) ...
            - sum(sum(ET*dth))/1000*(cellsize^2) + sum(sum(EX))/1000*(cellsize^2)...
            -  sum(sum(QsurR*dth))/1000*(cellsize^2)  -  sum(sum(QchR*dth))/1000*(cellsize^2) ...
            - sum(sum(Lk*dth))/1000*(cellsize^2)+ Qsub/1000*(cellsize^2)*dth ...
            - QoutS(i)/1000*(cellsize^2)  -QoutC(i)/1000*(cellsize^2) -Qout(i)/1000*(cellsize^2) ;% - QpointS(i)/1000*(cellsize^2) -Qpoint(i)/1000*(cellsize^2);  %%[m^3]
        CK2(i)=1000*CK2(i)/Area; %%[mm]
        
        
        Oave(i)= mean(reshape(O(O~=0),1,numel(O(O~=0))));
        Ostd(i)= std(reshape(O(O~=0),1,numel(O(O~=0))));
        
        S_ave(i)= (Oave(i)-Ohy(Yout,Xout))./(Osat(Yout,Xout)-Ohy(Yout,Xout));
        Upoint(i)= cellsize./tH_store(Yout,Xout); %%% Surface Velocity [m/s]
        UCh_point(i)= cellsize./t_store(Yout,Xout); %%% Surface Velocity [m/s]
        
    end
    
    % Plot soil moisture
    set(figure(1),'position',[100 100 500 800])
    subplot(211)
    imagesc((O-Ohy)./(Osat-Ohy))
    caxis([0 1])
    title(['Effective Saturation T = ' num2str(i) ' h'])
    colorbar
    axis equal 
    axis tight
    
    subplot(212)
    plot((2:Ndt),Qpoint(2:end)/1000*(cellsize^2)/dt_in_d,':m','Linewidth',1.5)
    hold on ; grid on;
    plot((2:Ndt),QpointS(2:end)/1000*(cellsize^2)/dt_in_d,':g','Linewidth',1.5)
    grid on ; hold on ;
    plot((2:Ndt),QpointC(2:end)/1000*(cellsize^2)/dt_in_d,'r','Linewidth',1.5)
    grid on ; hold on
    ylabel(['Q [m^3/s]'])
    ylim([0 20])
    pause(0.1)

end


figure(2)
subplot(211)
bar((1:Ndt),Prt)
ylabel('P [mm/h]')

subplot(212)
plot((2:Ndt),Qpoint(2:end)/1000*(cellsize^2)/dt_in_d,':m','Linewidth',1.5)
hold on ; grid on;
plot((2:Ndt),QpointS(2:end)/1000*(cellsize^2)/dt_in_d,':g','Linewidth',1.5)
grid on ; hold on ;
plot((2:Ndt),QpointC(2:end)/1000*(cellsize^2)/dt_in_d,'r','Linewidth',1.5)
grid on ; hold on
ylabel(['Q [m^3/s]'])