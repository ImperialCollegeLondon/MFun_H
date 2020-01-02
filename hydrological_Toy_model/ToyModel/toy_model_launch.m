%%%%%%%%%%%%%%%%%%%%%%%
%%%% PROVA ROUTING -- COMPLESSIVO --
%%%%%%%%%%%%%%%%%%%%%%%
%LASTN = maxNumCompThreads(1);
%%%  DTM    S    T      R     SN
%%%% DTM --> Digital Elevation Model
%%%% SN--> StreamNetwork
%%%  T ---> flow_matrix sparse
%%%  R --> flow direction D_inf
%%%% S--> fraction []
%load data_routing_mig
load DTM_V11; T = T_flow2; [R, S] = dem_flow(DTM,cellsize,cellsize);
%load data_routing2_biosphere2
%cellsize=10;  %%%% m
%xllcorner=10.5; yllcorner=10.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n]=size(DTM);
x=xllcorner:cellsize:(xllcorner+cellsize*(n-1));
y=yllcorner:cellsize:(yllcorner+cellsize*(m-1));
%%%%%%%%%%%%%%%%%%%%%%
MASK=ones(m,n); MASK(isnan(DTM))=0;
SLO = S; clear S %% fraction
SLO(SLO<0.0001)=0.0001;
SLO=SLO.*MASK;
SN(isnan(SN))=0;
%SLO(SN==1)=0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Matrice Zs [mm]
Zs=zeros(m,n);
Zs(:,:) = 1000.*MASK; %%  Soil Active depth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Osat =  (zeros(m,n) + 0.41).*MASK;
Ohy =   (zeros(m,n) + 0.057).*MASK;
Oel =   (zeros(m,n) + 0.11).*MASK;
%%%%%%%%%%%%%%%%%%
nvg= 2.28;
avg = -6.3; %% [1/m]
Ks=  (zeros(m,n) + 146).*MASK; %%% [mm/h]
mvg = (zeros(m,n) + (1 - 1/nvg)).*MASK;
Kbot=  (zeros(m,n) + 0.0).*MASK; %%% [mm/h]
%%%%%%%%%%%%%%%%%%%%%%%%%%
kF = 10000;%   %%% [h]  acquifer constant
mF= 100000; %320; %%% [mm]
f=(Osat-Ohy)./mF; %%% [1/mm]
%%%%%
OPT_UNSAT =0; %%% Option unsaturated 1 -- saturated bottom 0
CF=0; %% Acquifer Yes - No
aR=1; %%[-] anysotropy ratio
aT= 1000*(cellsize^2)/cellsize; %%[mm] Area/Contour length ratio
Area= (cellsize^2)*sum(sum(MASK)); %% [m^2]
Aacc = upslope_area(DTM, T)*cellsize^2; %%% [m^2]
%%%%%%
%%%%% dth --- >>>>
dt =360;%%%[s]
dth= dt/3600;% [h]
dti= 10; %%[s] Internal Time step for Surface Overland-flow Routing
dti2= 2; %%[s] Internal Time step for Surface Channel-flow Routing
%%%%%%%%%%%%%%%%%%%%%%%%
WC = cellsize*ones(m,n); %%% [m]  width channel
WC(SN==1)=0.0025*sqrt(Aacc(SN==1)); %% [m]
WC=WC.*MASK;
NMAN_C=SN*0.03; NMAN_H=0.10; NMAN_C=NMAN_C.*MASK;%%[s/(m^1/3)] manning coefficient
NMAN_H=NMAN_H.*MASK;
MRO = 0.0; %0.002; %%[m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vmax=zeros(m,n);
Vmax =Zs.*(Osat-Ohy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
Ndt=680;%7200;%8640;
%%%%%%%%%%%%%%%
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
Xout=8; Yout=1;
%Xout=255; Yout = 1;
%Xout=4; Yout=9;
%Xout =192; Yout =50;
SLO(Yout,Xout)=0.05;
%%%% Tracked point
%Xout =166; Yout =88;
Xout=8; Yout=18;
%%%%%%% Image plot
i_PR=98000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Precipitation
Prt=zeros(1,Ndt); %%% [mm/h]
Prt(2:91)=9.4;
%%% Evapotranspiration
Rv = 461.5 ;%% [J/kg K]
RH=0.6;
ETr=zeros(1,Ndt);
ETPt=zeros(1,Ndt);
ETPt(2:Ndt)=0.0;% [mm/h]  %%
%%% Initial Condition
QsurR=zeros(m,n);
QchR=zeros(m,n);
%%%%%%%%%%%%
Vf=zeros(1,Ndt);
Vf(1)=0;
%%%%%%% Intial Oi
%Oi=(Osat-0.15).*MASK;   O = Oi;
Oi=0.26.*MASK; O=Oi;
%Oi=unifrnd(0.08,0.12,m,n); O = Oi;
%Oi=unifrnd(0.06,0.10,m,n); O = Oi;
%%%%%%%%%%%%%%%%%%%
CK1=zeros(1,Ndt);
CK2=zeros(1,Ndt);

col='og';
%%%


tic ;
%%%%%%%%
bau = waitbar(0,'Time Step');
%%%%%%%%%%%%%%%%%%
for i=2:Ndt
    %%%%%%%%%%%%%%%%%%
    waitbar(i/Ndt,bau);
    %%% Matrix Flow
    Pr= Prt(i)*MASK + QsurR; %%% [mm/h]  Precipitation
    ET= ETPt(i)*MASK; %%% [mm/h]  evapotranspiration
    %BETA=ones(m,n);
    BETA=(O(:,:)-Ohy)./(Oel-Ohy);  BETA(BETA>1)=1; BETA(BETA<0)=0;
    BETA(isnan(BETA))=0; %%%%% evapotranspiration factor
    ET=ET.*BETA; %[mm/h]  evapotranspiration
    %%%%%%%%
    ETr(i) =  sum(sum(ET))*(cellsize^2/Area); %[mm/h]  evapotranspiration
    %%%%%
    Rh = Pr - Ks;  Rh(Rh<0)=0; %%% [mm/h] Horton Runoff
    If = Pr; If(If>Ks)=Ks(If>Ks); If=If*dth; %% [mm] Infiltration
    %%%%%%% First Layer %%%%
    VI= (O-Ohy).*Zs; %%% Volume First layer present [mm]
    V= VI + If;  %% Volume first update [mm]
    O = V./Zs + Ohy;  O(isnan(O))=0; %%
    O(O>Osat)=Osat(O>Osat); %% Soil Moisture first update []
    %%%%%%%%%%%%%
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
    %%%%%%%%%%%%%%%%%%%%%%%%
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
            (1-(1-((O-Ohy)./(Osat-Ohy)).^(1./mvg)).^mvg).^2; %%%%%% Unsaturated conductivity surface [mm/h]
        Lk(isnan(Lk))=0;
    else
        Lk=CF*Ks_b; %%%  [mm/h]
        Lk(isnan(Lk))=0;
    end
    %%%%%%%%%%%%%%
    Qi = kdtS.*To.*SLO/aT; %%% [mm/h] %%% Lateral Outflow
    V= V-Qi*dth - Lk*dth- ET*dth; %%% Volume second update [mm]
    %%%%%%%%
    Rd = V-Vmax; Rd(Rd<0)=0; % Dunne Runoff [mm]
    V=V-Rd; EX = -V.*(V<0); V(V<0)=0; %%% Volume Third update
    Rd=Rd/dth; %%% Dunne Runoff [mm/h]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Acquifer
    Vf(i) = Vf(i-1) + sum(sum(Lk*dth)) -sum(sum(EX)); %%% [mm] Totali
    if Vf(i)<0
        disp('Time Step too Long')
        break
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Qsur= Rd + Rh; %%% Surface flow [mm/h]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CK1(i)= sum(sum(VI))/1000*(cellsize^2)  -  sum(sum(V))/1000*(cellsize^2) + sum(sum(Pr*dth))/1000*(cellsize^2) ...
        - sum(sum(ET*dth))/1000*(cellsize^2) + sum(sum(EX))/1000*(cellsize^2) ...
        -  sum(sum(Qsur*dth))/1000*(cellsize^2) -  sum(sum(Qi*dth))/1000*(cellsize^2) - sum(sum(Lk*dth))/1000*(cellsize^2);  %%[m^3]
    CK1(i)=1000*CK1(i)/Area; %%[mm]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Subsurface Routing
    Qseep=Qi*dth.*(SN); %%[mm]
    Qi=Qi.*(1-SN); %% [mm/h]
    %%% SUBSURFACE VELOCITY
    [QiM]=Flow_Routing_Step2(DTM,T,Qi*dth); %%[mm]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    QoutS(i)= (sum(sum(Qi*dth))- sum(sum(QiM))); %%%  %%[mm]
    QpointS(i)=  Qi(Yout,Xout)*dth;%%%  %%[mm]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    %%%%%%%%%%%%%
    Qi_fall = 0*MASK;
    tH_store = 0*MASK;
    %%%%% OVERLAND FLOW
    Qsur = Qsur*dth; %%[mm]
    if sum(sum(Qsur))>0
        for jjj=1:1:dt/dti;
            %%%% SURFACE VELOCITY
            %WC= WC0 + 2*Y*cot(alpha); %%% [m]
            %t=(cellsize.*(WC.^0.4).*(NMAN.^0.6))./(((cellsize^2)*Qsur/3600/1000).^0.4.*(SLO).^0.3); %%% [s]
            YH = (Qsur*dth/1000);% .*cellsize./WC ; %% [m]
            %%%%% Ponding - Microroughness 
            YH = YH - MRO ; YH(YH<0)=0; 
            %%%%      
            tH= (cellsize.*NMAN_H)./(YH.^(2/3).*SLO.^0.5); %%% [s]
            tH(isnan(tH))=0;
            tH_store = tH_store + tH*dti/dt;
            kdt= dti./tH; %[]
            kdt(kdt>1)=1;
            %kdt(kdt>0)=1; %%% option no-routing
            kdt(isnan(kdt))=0;
            %%%% Surface Routing
            [QsurM]=Flow_Routing_Step2(DTM,T,kdt.*Qsur); %%[mm]
            QsurM2 = QsurM.*(1-SN); %%[mm]
            Qi_fall = Qi_fall + QsurM.*(SN); %% [mm]
            I_fall = sum(sum(QsurM.*(SN)));
            QsurM = QsurM2;
            QsurR = QsurM + (Qsur - kdt.*Qsur) ; %%[mm]
            %%%%%%%%%%
            Qout(i)= Qout(i) + (sum(sum(Qsur))- sum(sum(QsurR))-I_fall); %%%%% -- [mm]
            Qpoint(i)= Qpoint(i) + kdt(Yout,Xout)*Qsur(Yout,Xout); %%%%% -- [mm]
            %%%%
            Qsur=QsurR; %%[mm]
        end
    else
        YH =Qsur ; %% [m]
        tH= Inf*Qsur; %%% [s]
        Qout(i)= 0; %%%%% -- [mm]
        Qpoint(i)=0; %%%%% -- [mm]
        %%%%%%%%%%%%%%%%%%%%%
        QsurR=0; %%[mm]
        Qi_fall=0; %%[mm]
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    Qch = QchR*dth + Qseep + Qi_fall; %%% [mm]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t_store = 0*MASK;
    %%%
    if sum(sum(Qch))>0
        for jjj=1:1:dt/dti2;
            %%%%%% SURFACE VELOCITY %%%%%%%%%%%
            Y = (Qch/1000).*cellsize./WC ; %% [m]
            t= (cellsize.*NMAN_C)./(Y.^(2/3).*SLO.^0.5); %%% [s]
            t(isnan(t))=0;
            t_store = t_store + t*dti2/dt;
            kdt= dti2./t; %[]
            kdt(kdt>1)=1;
            kdt(isnan(kdt))=0;
            %%%% Surface Routing
            [QchM]=Flow_Routing_Step2(DTM,T,kdt.*Qch); %%[mm]
            QchR = QchM + (Qch - kdt.*Qch) ; %% [mm]
            %%%%%%%%%%%%%%
            QoutC(i)= QoutC(i) + (sum(sum(Qch))- sum(sum(QchR))); %%%%% -- [mm]
            QpointC(i)= QpointC(i) + kdt(Yout,Xout)*Qch(Yout,Xout); %%%%% -- [mm]
            %%%
            Qch=QchR; %%[mm]
        end
    else
        Y =Qch ; %% [m]
        t= Qch*Inf; %%% [s]
        QoutC(i)= 0; %%%%% -- [mm]
        QpointC(i)=0; %%%%% -- [mm]
        %%%%%%%%%%%%%%%%%%%%%
        QchR=0; %%[mm]
    end
    %%%%%%%%%%%%%%%%%%%
    QsurR = QsurR/dth; %%% [mm/h]
    QchR = QchR/dth; %%% [mm/h]
    QiM = QiM/dth; %%% [mm/h]
    %%%%%%%%%
    V = V + QiM*dth  ; %%% Volume fourth update [mm]
    Rd2 = V-Vmax; Rd2(Rd2<0)=0; % Dunne Runoff 2 [mm]
    V=V-Rd2; %%% Volume fifth update
    O= V./Zs + Ohy; O(isnan(O))=0; %% % Soil Moisture update []
    QsurR= QsurR + Rd2/dth;  %%% [mm/h]
    %%%%%%%%%%%%%%%%%%%%%
    %%%% Acquifer Routing
    Qsub= Vf(i)/kF; %%% [mm/h]
    Vf(i)=Vf(i)- Qsub*dth; %%[mm]
    QchR= QchR + SN*Qsub/sum(sum(SN));  %%% [mm/h]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CK2(i)= sum(sum(VI))/1000*(cellsize^2)  -  sum(sum(V))/1000*(cellsize^2) + sum(sum(Pr*dth))/1000*(cellsize^2) ...
        - sum(sum(ET*dth))/1000*(cellsize^2) + sum(sum(EX))/1000*(cellsize^2)...
        -  sum(sum(QsurR*dth))/1000*(cellsize^2)  -  sum(sum(QchR*dth))/1000*(cellsize^2) ...
        - sum(sum(Lk*dth))/1000*(cellsize^2)+ Qsub/1000*(cellsize^2)*dth ...
        - QoutS(i)/1000*(cellsize^2)  -QoutC(i)/1000*(cellsize^2) -Qout(i)/1000*(cellsize^2) ;% - QpointS(i)/1000*(cellsize^2) -Qpoint(i)/1000*(cellsize^2);  %%[m^3]
    CK2(i)=1000*CK2(i)/Area; %%[mm]
    %%%%%%%%%%%%%%%%%%%%
    %%%%%%
    %%%%%%%%%%%%%%%%%%%
    Oave(i)= mean(reshape(O(O~=0),1,numel(O(O~=0))));
    Ostd(i)= std(reshape(O(O~=0),1,numel(O(O~=0))));
    %Qi_std= std(reshape(QiM,1,numel(QiM)));
    %ET_std= std(reshape(ET,1,numel(ET)));
    %r_QO=corrcoef(reshape(O,1,numel(O)),reshape(QiM,1,numel(QiM)));
    %r_EO=corrcoef(reshape(O,1,numel(O)),reshape(ET,1,numel(ET)));
    %Cov_OQ(i)= r_QO(1,2)*Qi_std*Ostd(i);
    %Cov_OE(i)= r_EO(1,2)*ET_std*Ostd(i);
    X=reshape(O(O~=0),1,numel(O(O~=0)));
    edges=Ohy(Yout,Xout):(Osat(Yout,Xout)-Ohy(Yout,Xout))/50:Osat(Yout,Xout);
    pX = histc(X,edges)/length(X);
    pX=pX(pX~=0);
    ENT(i)=-sum(pX.*log(pX));
    ENT2(i)=(1/(3/4-1))*(1-sum(pX.^(3/4)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Helmholtz Free energy
    Se = (O-Ohy)./(Osat-Ohy);
    Psi = (1./avg).*((Se).^(-1./mvg)-1).^(1/nvg); %%% [m]
    FEN(i)=  9810*cellsize*cellsize*0.001*Zs(1,1)*(sum(sum(Psi.*O)) +  sum(sum((DTM-0.001*Zs/2-min(min(DTM))).*O))); %%[J]
    %%%
    EPE(i)= -Rv*ETr(i)/3.6*log(RH); % Entropy production ET [mW/m^2 K]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S_ave(i)= (Oave(i)-Ohy(Yout,Xout))./(Osat(Yout,Xout)-Ohy(Yout,Xout));
    Upoint(i)= cellsize./tH_store(Yout,Xout); %%% Surface Velocity [m/s]
    UCh_point(i)= cellsize./t_store(Yout,Xout); %%% Surface Velocity [m/s]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%
    if i== i_PR || i == Ndt
        figure(10000+i)
        imagesc(x,y,QsurR*dth); axis xy  ; axis equal
        title('Overland flow [mm] after Ndt time step'); colorbar;
        figure(20000+i)
        imagesc(x,y,QchR*dth); axis xy  ; axis equal
        title('Discharge [mm] after Ndt time step'); colorbar;
        figure(30000+i)
        imagesc(x,y,cellsize./tH_store); axis xy  ; axis equal
        title('Surface flow velocity hillslope [m/s] after Ndt time step'); colorbar;
        figure(40000+i)
        imagesc(x,y,SN*cellsize./t_store); axis xy  ; axis equal
        title('Surface flow velocity channel [m/s] after Ndt time step'); colorbar;
        figure(50000+i)
        imagesc(x,y,(O-Ohy)./(Osat-Ohy)); axis xy  ; axis equal
        title('Effective Saturation'); colorbar;
        figure(60000+i)
        imagesc(x,y,O); axis xy  ; axis equal
        title('Soil Moisture'); colorbar; %caxis([Ohy(1,1) Osat(1,1)]); %caxis([0 max(max(O(:,:,2)))]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
close(bau);
TIME_elapsed = toc/Ndt*1000 %%% [ms /step ]
disp('[ms/step]')


%%%% Check Mass
VI= (Oi-Ohy).*Zs; %%% Volu
CK=   sum(sum(VI))/1000*(cellsize^2)  -  sum(sum(V))/1000*(cellsize^2) ...
    + sum(sum(Prt(2:Ndt)*dth))/1000*Area - sum(sum(ETr(2:Ndt)*dth))/1000*Area ...
    -  sum(sum(QsurR*dth))/1000*(cellsize^2)   -  sum(sum(QchR*dth))/1000*(cellsize^2)  + ...
    (Vf(1)-Vf(end))/1000*(cellsize^2)-sum(QoutS)/1000*(cellsize^2)- sum(QoutC)/1000*(cellsize^2) - sum(Qout)/1000*(cellsize^2);  %%% [m^3]


%%%%%%%%%%%%%%%%%%%%% FINAL FIGURE
figure(4)
plot((2:Ndt)*dth,Qpoint(2:end)/1000*(cellsize^2)/dt,':m','Linewidth',1.5)
hold on ; grid on;
plot((2:Ndt)*dth,QpointS(2:end)/1000*(cellsize^2)/dt,':g','Linewidth',1.5)
grid on ; hold on ;
plot((2:Ndt)*dth,QpointC(2:end)/1000*(cellsize^2)/dt,'r','Linewidth',1.5)
grid on ; hold on ;
xlabel('Time [h]'); ylabel('Discharge [m^3/s]');
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
subplot(2,1,1)
plot((2:Ndt)*dth,Qout(2:end)/dth*(cellsize^2)/Area,':m','Linewidth',1.5)
hold on ; grid on;
plot((2:Ndt)*dth,(QoutS(2:end))/dth*(cellsize^2)/Area,':g','Linewidth',1.5)
grid on ; hold on
plot((2:Ndt)*dth,(QoutC(2:end))/dth*(cellsize^2)/Area,'k','Linewidth',1.5)
xlabel('Time [h]'); ylabel('Discharge [mm/h]');
subplot(2,1,2)
plot((2:Ndt)*dth,Qout(2:end)/1000*(cellsize^2)/dt,':m','Linewidth',1.5)
hold on ; grid on;
plot((2:Ndt)*dth,(QoutS(2:end))/1000*(cellsize^2)/dt,':g','Linewidth',1.5)
grid on ; hold on ;
plot((2:Ndt)*dth,(QoutC(2:end))/1000*(cellsize^2)/dt,'k','Linewidth',1.5)
xlabel('Time [h]'); ylabel('Discharge [m^3/s]');




figure(6)
plot((2:Ndt)*dth,Vf(2:end)/1000*(cellsize^2)/kF/3600,'r','Linewidth',3)
hold on ; grid on;
xlabel('Time [h]'); ylabel('Acquifer Flow [m^3/s]');
figure(7)
plot((2:Ndt)*dth,S_ave(2:end),'r','Linewidth',3)
hold on ; grid on;
xlabel('Time [h]'); ylabel('Effective Saturation Basin');
figure(8)
plot((2:Ndt)*dth,Upoint(2:end),'--g','Linewidth',3)
hold on ; grid on;
plot((2:Ndt)*dth,UCh_point(2:end),'r','Linewidth',3)
xlabel('Time [h]'); ylabel('Surface Velocity');
figure(80)
plot3((2:Ndt)*dth,Oave(2:end),Ostd(2:Ndt)./Oave(2:end),'ob','Linewidth',2)
hold on ; grid on;
ylabel('Average Soil Moisture'); zlabel('Cv -- Soil Moisture');
xlabel('Time [h]');


%%%%%%%%%%%%%%%%%%%%%%%%%

figure(82)
plot(Oave(2:end),Ostd(2:Ndt)./Oave(2:end),col,'Linewidth',2)
hold on ; grid on;
xlabel('Average Soil Moisture'); ylabel('Cv -- Soil Moisture');
%%%%%%%%%%%%%
%   figure(3)
%   plot((2:Ndt)*dth,Ostd(2:Ndt),col(2),'Linewidth',2)
%   hold on ; grid on;
%   plot((2:Ndt)*dth,Oave(2:end),strcat('--',col(2)),'Linewidth',2)
%   xlabel('Time [h]');ylabel('Ave/Std -- Soil Moisture');
%   figure(4)
%   plot((2:Ndt)*dth,Cov_OQ(2:end),col(2),'Linewidth',2)
%   hold on ; grid on;
%   plot((2:Ndt)*dth,Cov_OE(2:end),strcat('--',col(2)),'Linewidth',2)
%   xlabel('Time [h]');ylabel('Cov \theta-Qi  \theta-ET ');

figure(85)
plot((2:Ndt)*dth,ENT(2:end),col(2),'Linewidth',3)
hold on ; grid on;
xlabel('Time [h]'); ylabel('ENTROPY OF THE SYSTEM');
figure(86)
subplot(3,1,1)
plot((2:Ndt)*dth,FEN(2:end),col(2),'Linewidth',2)
hold on ; grid on;
xlabel('Time [h]'); ylabel('[J]') ; title('Helmholtz Free Energy');
subplot(3,1,2)
plot((122:Ndt)*dth,-1000*diff(FEN(121:end))/60/Area/288,col(2),'Linewidth',2)
hold on ; grid on;
xlabel('Time [h]'); ylabel('[mW/m^2 K]') ; title('Entropy production by redistribution ');
subplot(3,1,3)
plot((2:Ndt)*dth,EPE(2:end),col(2),'Linewidth',2)
hold on ; grid on;
xlabel('Time [h]'); ylabel('[mW/m^2 K]') ; title('ET Entropy production');

figure(87)
plot(Oave(2:end),FEN(2:end),col,'Linewidth',1.5)
hold on ; grid on
xlabel('Average Soil Moisture'); ylabel('[J]') ; title('Helmholtz Free Energy');



ANSW_U = 0;
if ANSW_U == 1
    figure(100000+i)
    subplot(1,3,1)
    surf(x,y,DTM,QchR*dth); axis xy  ; axis equal
    title('Discharge [mm]'); colorbar;
    xlabel('[m]');ylabel('[m]');zlabel('[m]')
    subplot(1,3,2)
    surf(x,y,DTM,cellsize./t); axis xy  ; axis equal
    title('Surface flow velocity [m/s]'); colorbar;
    xlabel('[m]');ylabel('[m]');zlabel('[m]')
    subplot(1,3,3)
    surf(x,y,DTM,(O-Ohy)./(Osat-Ohy)); axis xy  ; axis equal
    title('Effective Saturation [-]'); colorbar;
    xlabel('[m]');ylabel('[m]');zlabel('[m]')
end



