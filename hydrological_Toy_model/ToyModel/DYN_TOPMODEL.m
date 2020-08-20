

function [OUTPUT] = DYN_TOPMODEL(TOPO_DAT,HYDR_DATA,METEO_DAT,SIM_PARAM,OUT_PARAM)


%% INPUT PREPROCESSING
DTM = TOPO_DAT.DTM;
T = TOPO_DAT.T;
SN = TOPO_DAT.SN;
WC = TOPO_DAT.WC;
MASK = TOPO_DAT.MASK;
Zs = TOPO_DAT.Zs; 
cellsize = TOPO_DAT.cellsize;

Osat = HYDR_DATA.Osat;
Ohy = HYDR_DATA.Ohy;
Oel = HYDR_DATA.Oel;
Ks = HYDR_DATA.Ks;
aR = HYDR_DATA.aR; %[-] anisotropy ratio
mvg = HYDR_DATA.mvg;
Kbot = HYDR_DATA.Kbot;
kF = HYDR_DATA.kF;%  [h]  acquifer constant
mF= HYDR_DATA.mF;% [mm]
NMAN_C=HYDR_DATA.NMAN_C; 
NMAN_H=HYDR_DATA.NMAN_H; 
MRO = HYDR_DATA.MRO; %0.002; %%[m]

f=(Osat-Ohy)./mF;

Prt = METEO_DAT.Prt;
ETPt = METEO_DAT.ETPt;
D_h = METEO_DAT.D_h;

dt = SIM_PARAM.dt;
dti = SIM_PARAM.dti;
dti2 = SIM_PARAM.dti2;
OPT_UNSAT = SIM_PARAM.OPT_UNSAT;
CF = SIM_PARAM.CF; % Acquifer Yes - No
Oi = SIM_PARAM.Oi;
pl = SIM_PARAM.pl;
    
Xout = OUT_PARAM.Xout;
Yout = OUT_PARAM.Yout;
T_map_s = OUT_PARAM.T_map_s;

%% Simulation Initialization
dth= dt/3600;% [h]
[~, S] = dem_flow(DTM,cellsize,cellsize);
SLO = S; 
SLO(SLO<0.0001)=0.0001;
SLO=SLO.*MASK;

np = numel(Xout);

NMAN_C=NMAN_C.*MASK;%[s/(m^1/3)] manning coefficient
NMAN_H=NMAN_H.*MASK;

aT= 1000*(cellsize^2)/cellsize; %[mm] Area/Contour length ratio
Area= (cellsize^2)*sum(sum(MASK)); % [m^2]

Ndt = length(Prt);
[m,n] = size(DTM);

Vmax = zeros(m,n);
Vmax =Zs.*(Osat-Ohy);
Qout=zeros(1,Ndt);
QoutS=zeros(1,Ndt);
Qpoint=zeros(Ndt,np);
QpointS=zeros(Ndt,np);
QoutC=zeros(1,Ndt);
QpointC=zeros(Ndt,np);
Upoint=zeros(Ndt,np);
UCh_point=zeros(1,Ndt);
S_ave=zeros(1,Ndt);
Oave=zeros(1,Ndt);
Ostd=zeros(1,Ndt);
QsurR=zeros(m,n);
QchR=zeros(m,n);
Vf=zeros(1,Ndt);
Vf(1)=0;
O=Oi;
CK1=zeros(1,Ndt);
CK2=zeros(1,Ndt);
OUTPUT.O = cell(1,ceil(Ndt/T_map_s));
t_store = zeros(m,n);
s_c = 1;

DTM_nml = numel(DTM);
[mm,nn]=size(DTM);
XX = -T + speye(mm*nn,nn*mm);

SNl = logical(SN);
Y = zeros(size(DTM));
tH = 0*DTM;
t = 0*DTM;
            
CHR_m = (cellsize.*NMAN_C)./(SLO.^0.5);
SR_m = (cellsize.*NMAN_H)./(SLO.^0.5);
CHR_m2 = cellsize./WC/1000;
    
for i=2:Ndt
    
    disp(i)
    
    % Matrix Flow
    Pr= Prt(i)*MASK + QsurR; % [mm/h]  Precipitation
    ET= ETPt(i)*MASK; % [mm/h]  evapotranspiration
    %BETA=ones(m,n);
    BETA=(O(:,:)-Ohy)./(Oel-Ohy);  BETA(BETA>1)=1; BETA(BETA<0)=0;
    BETA(isnan(BETA))=0; % evapotranspiration factor
    ET=ET.*BETA; %[mm/h]  evapotranspiration
    
    ETr(i) =  sum(sum(ET))*(cellsize^2/Area); %[mm/h]  evapotranspiration
    Rh = Pr - Ks;  Rh(Rh<0)=0; % [mm/h] Horton Runoff
    
    If = Pr;
    If(If>Ks)=Ks(If>Ks);
    If=If*dth; % [mm] Infiltration
    
    % First Layer
    VI= (O-Ohy).*Zs; % Volume First layer present [mm]
    V= VI + If;  % Volume first update [mm]
    O = V./Zs + Ohy;  O(isnan(O))=0;
    O(O>Osat)=Osat(O>Osat); % Soil Moisture first update []
    
    Ko = Ks.*(((O-Ohy)./(Osat-Ohy)).^(1/2)).*...
        (1-(1-((O-Ohy)./(Osat-Ohy)).^(1./mvg)).^mvg).^2; %%%%%% Unsaturated conductivity surface [mm/h]
    Ko(isnan(Ko))=0; 
    if OPT_UNSAT == 1
        tS = (Osat).*1000*cellsize./Ko; % [h]
        tS(isnan(tS))=0; 
    else
        tS = (Osat).*1000*cellsize./Ks; % [h]
        tS(isnan(tS))=0; 
    end
    
    kdtS=1;
    
    if OPT_UNSAT == 1
        To = (aR.*Ko./f).*(1 -exp(-f.*Zs)); % Trasmsissivity [mm^2/h]
        To(isnan(To))=0;
    elseif OPT_UNSAT == 2
        Z= V./(Osat-Ohy);% [mm]
        Z(isnan(Z))=0;
        To= Ks.*(Z); % [mm^2/h]
    elseif OPT_UNSAT == 3
        % Original Topmodel formulation
        
    end
    
    Ks_Z = Ks.*exp(-f.*Zs); Ks_Z(Ks_Z<Ks/10)= Ks(Ks_Z<Ks/10)/10; % [mm/h] Saturated Conductivity Bottom
    Ks_b =  exp((log(Kbot) + log(Ks_Z))/2);%  [mm/h] Conductivity between layer n and bedrock
    if OPT_UNSAT == 1
        Lk= CF*Ks_b.*(((O-Ohy)./(Osat-Ohy)).^(1/2)).*...
            (1-(1-((O-Ohy)./(Osat-Ohy)).^(1./mvg)).^mvg).^2; % Unsaturated conductivity surface [mm/h]
        Lk(isnan(Lk))=0;
    else
        Lk=CF*Ks_b; %  [mm/h]
        Lk(isnan(Lk))=0;
    end
    
    Qi = kdtS.*To.*SLO/aT; % [mm/h] %%% Lateral Outflow
    V= V-Qi*dth - Lk*dth- ET*dth; % Volume second update [mm]
    
    Rd = V-Vmax; Rd(Rd<0)=0; % Dunne Runoff [mm]
    V=V-Rd; EX = -V.*(V<0); V(V<0)=0; % Volume Third update
    Rd=Rd/dth; % Dunne Runoff [mm/h]
    
    % Acquifer
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
    
    % Subsurface Routing
    Qseep=Qi*dth.*(SN); %%[mm]
    Qi=Qi.*(1-SN); %% [mm/h]
    
    Qsurl = reshape(Qi*dth,DTM_nml,1);
    QiM =(XX)*Qsurl;
    QiM  = reshape(QiM , [mm nn]);
    
    QoutS(i)= (sum(sum(Qi*dth))- sum(sum(QiM))); %%%  %%[mm]
    
    for p_i = 1:np
        QpointS(i,p_i)=  Qi(Yout(p_i),Xout(p_i))*dth;%%%  %%[mm]
    end
    
    Qi_fall = 0*MASK;
    tH_store = 0*MASK;
    % OVERLAND FLOW
    Qsur = Qsur*dth; %%[mm]
    Qsur = sparse(Qsur);
    
    if any(Qsur(:))
        for jjj=1:1:dt/dti
            % SURFACE VELOCITY
            YH = (Qsur*dth/1000);% .*cellsize./WC ; %% [m]
            % Ponding - Microroughness
            
            YH = YH - MRO ; YH(YH<0)=0;
%             tH= (cellsize.*NMAN_H)./(YH.^(2/3).*SLO.^0.5); %%% [s]
            tH= SR_m./(YH.^(2/3)); %%% [s]

            tH(isnan(tH))=0;
            tH_store = tH_store + tH*dti/dt;
            kdt= dti./tH; %[]
            kdt(kdt>1)=1;
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
            for p_i = 1:np
                Qpoint(i,p_i)= Qpoint(i,p_i) + kdt(Yout(p_i),Xout(p_i))*Qsur(Yout(p_i),Xout(p_i)); %%%%% -- [mm]
            end
            
            Qsur=QsurR; %[mm]
            Qsur = sparse(Qsur);
        end
    else
        YH =Qsur ; % [m]
        tH= Inf*Qsur; % [s]
        Qout(i)= 0; % -- [mm]
        Qpoint(i)=0; %-- [mm]
        
        QsurR=0; %[mm]
        Qi_fall=0; %[mm]
    end
    
    Qch = QchR*dth + Qseep + Qi_fall; % [mm]
%     Qch = sparse(Qch);
    t_store(:) = 0;

    if any(Qch(:))
        for jjj=1:1:dt/dti2
            
            % SURFACE VELOCITY 
            Y(SNl) = (Qch(SNl)).*CHR_m2(SNl) ; %% [m]
%             t(SNl)= (cellsize.*NMAN_C(SNl))./(Y(SNl).^(2/3).*SLO(SNl).^0.5); %%% [s]
            t(SNl) = CHR_m(SNl)./Y(SNl).^(2/3);
            t(isnan(t))=0;
            t_store = t_store + t*dti2/dt;
            kdt= dti2./t; %[]
            kdt(kdt>1)=1;
            kdt(isnan(kdt))=0;
            
            Qchl = reshape(kdt.*Qch,DTM_nml,1);
            QchM=(XX)*(Qchl);
            QchM = reshape(QchM, [mm nn]);
            QchR = QchM + (Qch - kdt.*Qch) ; % [mm]
            
            QoutC(i)= QoutC(i) + (sum(sum(Qch))- sum(sum(QchR))); % -- [mm]
            for p_i = 1:np
                QpointC(i,p_i)= QpointC(i,p_i) + kdt(Yout(p_i),Xout(p_i))*Qch(Yout(p_i),Xout(p_i)); %%%%% -- [mm]
            end
            Qch=QchR; %[mm]
        end
    else
        Y =Qch ; % [m]
        t= Qch*Inf; % [s]
        QoutC(i)= 0; % -- [mm]
        QpointC(i)=0; % -- [mm]
        QchR=0; %[mm]
    end
    
    %% addaptive time step addition
    
    t1 = nanmin(nanmin(tH(SN==0)));
    t2 = nanmin(nanmin(t(SN==1)));
    
    dti  = nanmax(t1/5,30);
    dti2 = nanmax(t2/5,2);
    
    dti = dt/ceil(dt/dti);
    dti2 = dt/ceil(dt/dti2);
    %%
    
    QsurR = QsurR/dth; % [mm/h]
    QchR = QchR/dth; % [mm/h]
    QiM = QiM/dth; % [mm/h]
    
    V = V + QiM*dth  ; % Volume fourth update [mm]
    Rd2 = V-Vmax; Rd2(Rd2<0)=0; % Dunne Runoff 2 [mm]
    V=V-Rd2; % Volume fifth update
    O= V./Zs + Ohy; O(isnan(O))=0; % Soil Moisture update []
    QsurR= QsurR + Rd2/dth;  % [mm/h]
    
    % Acquifer Routing
    Qsub= Vf(i)/kF; % [mm/h]
    Vf(i)=Vf(i)- Qsub*dth; %[mm]
    QchR= QchR + SN*Qsub/sum(sum(SN));  % [mm/h]
    
    CK2(i)= sum(sum(VI))/1000*(cellsize^2)  -  sum(sum(V))/1000*(cellsize^2) + sum(sum(Pr*dth))/1000*(cellsize^2) ...
        - sum(sum(ET*dth))/1000*(cellsize^2) + sum(sum(EX))/1000*(cellsize^2)...
        -  sum(sum(QsurR*dth))/1000*(cellsize^2)  -  sum(sum(QchR*dth))/1000*(cellsize^2) ...
        - sum(sum(Lk*dth))/1000*(cellsize^2)+ Qsub/1000*(cellsize^2)*dth ...
        - QoutS(i)/1000*(cellsize^2)  -QoutC(i)/1000*(cellsize^2) -Qout(i)/1000*(cellsize^2) ;% - QpointS(i)/1000*(cellsize^2) -Qpoint(i)/1000*(cellsize^2);  %%[m^3]
    CK2(i)=1000*CK2(i)/Area; %[mm]
    
    
    Oave(i)= mean(reshape(O(O~=0),1,numel(O(O~=0))));
    Ostd(i)= std(reshape(O(O~=0),1,numel(O(O~=0))));
    
    for p_i = 1:np
        S_ave(i,p_i)= (Oave(i)-Ohy(Yout(p_i),Xout(p_i)))./(Osat(Yout(p_i),Xout(p_i))-Ohy(Yout(p_i),Xout(p_i)));
        Upoint(i,p_i)= cellsize./tH_store(Yout(p_i),Xout(p_i)); % Surface Velocity [m/s]
        UCh_point(i,p_i)= cellsize./t_store(Yout(p_i),Xout(p_i)); % Surface Velocity [m/s]
    end
    
    if mod(i,T_map_s) == 0
        OUTPUT.O{s_c} = O;
        OUTPUT.T_o = D_h(i);
        s_c = s_c + 1;
        
    end
    
    if pl == 1
        
        set(figure(1),'position',[100 100 500 1000])
        subplot(211)
        imagesc((O-Ohy)./(Osat-Ohy))
        caxis([0 1])
        title(['Effective Saturation T = ' num2str(i) ' h'])
        colorbar
        axis equal
        axis tight
        
        subplot(212)
        plot((2:Ndt),Qpoint(2:end)/1000*(cellsize^2)/dt,':m','Linewidth',1.5)
        hold on ; grid on;
        plot((2:Ndt),QpointS(2:end)/1000*(cellsize^2)/dt,':g','Linewidth',1.5)
        grid on ; hold on ;
        plot((2:Ndt),QpointC(2:end)/1000*(cellsize^2)/dt,'r','Linewidth',1.5)
        grid on ; hold on
        ylabel(['Q [m^3/s]'])
        ylim([0 20])
        pause(0.1)
        
    end
end

OUTPUT.QpointC = QpointC;
OUTPUT.Qpoint = Qpoint;
OUTPUT.QpointS = QpointS;

end

