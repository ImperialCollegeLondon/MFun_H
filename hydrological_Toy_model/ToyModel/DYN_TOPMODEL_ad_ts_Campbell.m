function [OUTPUT] = DYN_TOPMODEL_ad_ts_Campbell(TOPO_DAT,HYDR_DATA,METEO_DAT,SIM_PARAM,OUT_PARAM)
% DYN_TOPMODEL_ad_ts provide a simplified model for rr modelling.
%
% Input Matrix type:
%  ETPt: Matrix<double>
%
%
% NOTE:
%  Campbell/Burdine formula were used for calculating Hydr.Conductivity &
%  Retention Curve;


% OUTPUT:
%   OUTPUT.O:      Soil Moisture
%   OUTPUT.T_o:
%   OUTPUT.QpointC:Channel Flow
%   OUTPUT.Qpoint: Overland Flow
%   OUTPUT.QpointS:Subsurface Flow


% Change all input time series into spatial version;
% Play with kF, anisotropy ratio
% Use data from Soil Grids;
% https://soilgrids.org/#!/?layer=ORCDRC_M_sl2_250m&vector=1
try
tag_Animation = 0;

%% INPUT PREPROCESSING
DTM = TOPO_DAT.DTM;
T = TOPO_DAT.T;
SN = TOPO_DAT.SN;
WC = TOPO_DAT.WC;
MASK = TOPO_DAT.MASK;
Zs = TOPO_DAT.Zs;
cellsize = TOPO_DAT.cellsize;

Osat = HYDR_DATA.Osat;
Osat(Osat == 0) = 1;
Ohy = HYDR_DATA.Ohy;
Oel = HYDR_DATA.Oel;
Ks = HYDR_DATA.Ks;
aR = HYDR_DATA.aR; %[-] anisotropy ratio
L = HYDR_DATA.L;% lambda
mvg = HYDR_DATA.mvg;

Kbed = HYDR_DATA.Kbed;
kF = HYDR_DATA.kF;%  [h]  acquifer constant 
bgw = HYDR_DATA.bgw;

Vmax = Zs.*(Osat-Ohy);

mF= HYDR_DATA.mF;% [mm]% calibrated using the observed data;
NMAN_C=HYDR_DATA.NMAN_C;
NMAN_H=HYDR_DATA.NMAN_H;
MRO = HYDR_DATA.MRO; %0.002; %%[m]


try 
   Pe = HYDR_DATA.Psi_ae; 
catch
    
end

f = (Osat-Ohy)./mF;

Prt = METEO_DAT.Prt;
PET = METEO_DAT.PET;
D_h = METEO_DAT.D_h;

dt = SIM_PARAM.dt;
dti = SIM_PARAM.dti;
dti2 = SIM_PARAM.dti2;
OPT_UNSAT = SIM_PARAM.OPT_UNSAT;
CF = SIM_PARAM.CF; % Acquifer Yes - No
Oi = SIM_PARAM.Oi;
pl = SIM_PARAM.pl;
GA = SIM_PARAM.GA;
try
    vf_ini = SIM_PARAM.vf_ini;
catch
    vf_ini = 2700;
end


Xout = OUT_PARAM.Xout;
Yout = OUT_PARAM.Yout;
T_map_s = OUT_PARAM.T_map_s;

%% Simulation Initialization
dth= dt/3600;% [h]
[~, S] = dem_flow(DTM,cellsize,cellsize);
SLO = S;
SLO(SLO<0.0001) = 0.0001;
SLO = SLO.*MASK;

np = numel(Xout);

NMAN_C = NMAN_C.*MASK;%[s/(m^1/3)] manning coefficient
NMAN_H = NMAN_H.*MASK;

aT = 1000*(cellsize^2)/cellsize; %[mm] Area/Contour length ratio
Area = (cellsize^2)*sum(sum(MASK)); % [m^2]

Ndt = Prt.length; % by Yuting

[m,n] = size(DTM);

[Qout,QoutS,QoutC,UCh_point,S_ave,Oave,Ostd,Vf] = deal(zeros(1,Ndt));
[Qpoint,QpointS,QpointC,Upoint] = deal(zeros(Ndt,np));
[QsurR,QchR]=deal(zeros(m,n));

if CF == 1
    Vf(1) = vf_ini;
else
    Vf(1) = 0;
end


O = Oi;
CK1 = zeros(1,Ndt);
CK2 = zeros(1,Ndt);
ETr = zeros(1,Ndt);

OUTPUT.O = cell(1,ceil(Ndt/T_map_s));
OUTPUT.SatSoil = NaN(Ndt,1);
t_store = zeros(m,n);
s_c = 1;

DTM_nml = numel(DTM);
[mm,nn] = size(DTM);
XX = -T + speye(mm*nn,nn*mm);

SNl = logical(SN);
Y = zeros(size(DTM));
tH = 0*DTM;
t = 0*DTM;

CHR_m = (cellsize.*NMAN_C)./(SLO.^0.5);
SR_m = (cellsize.*NMAN_H)./(SLO.^0.5);
CHR_m2 = cellsize./WC/1000;


if size(PET.ETPt) ~= 1
    fprintf('PET with %d h resolution is fully distributed\n',PET.scale);
end

[gwDisPart,leakage,ex] = deal(NaN(Ndt,1));

if size(PET.ETPt) == 1
    fprintf('PET: spatially uniform.\n');
else
    fprintf('PET: spatially distributed.\n');
end
if numel(size(Prt.rain)) == 3
    fprintf('RAINFALL: spatially distributed.\n');
else
    fprintf('RAINFALL: spatially uniform.\n');
end

% Green & Ampt related.
dryTag = Ohy*4;% # saving consecutive dry hours for each grid;
F = Ohy*0.0001; % the cumulative amount of water that has infiltrated;
If = Ohy*0; % infiltration rate computed from Green & Ampt;
psi_f = Ohy*0.0001;% matrix pressure at the wetting front;
O_ini = Oi;% iniitial moisture content before infiltration began;
KS_e = Ks;
Fp = Ohy*0;% threshold for Green & Ampt equation;


% tag: each grid is cat or not;

catTag = MASK;
catTag(catTag == 0) = NaN;


for i = 2:Ndt
    
    % disp(i);
    
    % INPUT: SPATIALLY DISTRIBUTED METEO + SOIL DATA for each step. by Yuting
    % [mm/h]  Precipitation <double>
    if numel(size(Prt.rain)) == 3
        Pr = squeeze(Prt.rain(:,:,i)).*MASK; 
        % Pr = RECRAIN(i).*MASK + QsurR; 
    else
        Pr = Prt.rain(i).*MASK; % [mm/h]  Precipitation<Matrix>
    end
    Pr = double(Pr);
    
    if GA == 1
        try
            
            DRYTHR = 3;% 3 % making it bigger can create more peaks
            
            % UPDATE dry hours Matrix
            % single point version:
            dryI = Pr < 0.002/24; % threshold for dry hours
            
            % spatial version:
            % coverRatio = sum(MASK(:) == 1 & Pr(:) > 0)/sum(MASK(:) == 1);
            % dryI = (coverRatio+Pr*0)<0.15;
            
            % case 1: COMPUTE psi_f, Oi for starting event
            newE = dryTag > DRYTHR & dryI == 0;
            % UPDATE Psi-f
            O33 = Oel;
            gw = 9810; %% specific weight water [N/m^3]
            B = 1./L;
            A = exp(log(33)+B.*log(O33)); % Coefficient of moisture tension
            if O < O33
                Psi_aux = A.*(O.^-B); %% [kPa]
            else
                Psi_aux = 33-((O-O33).*(33-Pe)./(Osat-O33));%% [kPa]
            end
            psi_f(newE) = 1000*1000*Psi_aux(newE)/(gw); %%[mm]% Tension at O
            O_ini(newE) = O(newE);
            F(newE) = 0.0000001;
            
            % case 2: UPDATE F for started event
            conE = ~(dryTag > DRYTHR & dryI == 0);
            F(conE) = F(conE) + If(conE);
            
            % UPDATE dry Hours #
            dryTag(dryI) = dryTag(dryI)+1;
            dryTag(~dryI) = 0;
            
        catch
            1;
        end
    end
    
    Pr = Pr + QsurR;
    
    % [mm/h]  evapotranspiration <single> -> <double>
    if size(PET.ETPt) == 1
        ET = double(squeeze(PET.ETPt(ceil(i/PET.scale))))/PET.scale.*MASK;
    else
        ET = double(squeeze(PET.ETPt(:,:,ceil(i/PET.scale))))/PET.scale.*MASK; 
    end
    
    % BETA: refer to [Bergström, 1992]
    % BETA = ones(m,n);
    BETA = (O(:,:)-Ohy)./(Oel-Ohy);  BETA(BETA>1) = 1;  BETA(BETA<0) = 0;
    BETA(isnan(BETA)) = 0; % soil water availability factor
    ET = ET.*BETA; %[mm/h]  evapotranspiration<Matrix>
    
    ETr(i) = sum(sum(ET))*(cellsize^2/Area); %[mm/h]  evapotranspiration
    Rh = Pr - Ks;  
    Rh(Rh<0) = 0; % [mm/h] Horton Runoff
    
    if GA == 1
        
        %%% Psi_ae [mm]
        % KS_e = Ks.*((Psi_ae.*(1- (O-Ohy)./(Osat-Ohy)).*Osat)./O./Zs +1);
        
        % Green-Ampt Function (\Cornell 2011\ http://www.hydrology.bee.cornell.edu/BEE3710Handouts/GreenAmpt.pdf)
        KS_e = Ks.*(1 + abs(psi_f).*(Osat-O_ini)./F);
        
        If = Pr;
        If(If>KS_e) = KS_e(If>KS_e);
        If = If*dth; % [mm] Infiltration
        
    elseif GA == 0
        
        If = Pr;
        If(If>Ks) = Ks(If>Ks);
        If = If*dth; % [mm] Infiltration
        
    elseif GA == 2
        
        
        
    end
    
    % First Layer
    VI = (O-Ohy).*Zs; % Volume First layer present [mm]
    V = VI + If;  % Volume first update [mm]
    O = V./Zs + Ohy;  O(isnan(O))=0;
    O(O > Osat) = Osat(O > Osat); % Soil Moisture first update []
    
    Ko = Ks.*(O./Osat).^(3+2./L); % Brooks & Corey Unsaturated conductivity surface [mm/h]
    %     Ko = Ks.*(((O-Ohy)./(Osat-Ohy)).^(1/2)).*...
    %         (1-(1-((O-Ohy)./(Osat-Ohy)).^(1./mvg)).^mvg).^2;
    if sum(isnan(Ko)) > 1
    disp('111');
    end
    
    Ko(isnan(Ko))=0;
    if OPT_UNSAT == 1
        tS = (Osat).*1000*cellsize./Ko; % [h]
        tS(isnan(tS)) = 0;
    else
        tS = (Osat).*1000*cellsize./Ks; % [h]
        tS(isnan(tS)) = 0;
    end
    
    kdtS=1;
    
    if OPT_UNSAT == 1
        To = (aR.*Ko./f).*(1 -exp(-f.*Zs)); % Trasmsissivity [mm^2/h]
        To(isnan(To)) = 0;
    elseif OPT_UNSAT == 2
        Z = V./(Osat-Ohy);% [mm]
        Z(isnan(Z)) = 0;
        To = Ks.*(Z); % [mm^2/h]
    elseif OPT_UNSAT == 3
        % Original Topmodel formulation
        %         ...........
    end
    
    Ks_Z = Ks.*exp(-f.*Zs); % [mm/h] Saturated Conductivity Bottom
    Ks_Z(Ks_Z<Ks/10)= Ks(Ks_Z<Ks/10)/10;

    Ks_b = exp((log(Kbed) + log(Ks_Z))/2);% [mm/h] Conductivity between layer n and bedrock
    
    if OPT_UNSAT == 1
        Lk = CF*Ks_b.*(O./Osat).^(3+2./L);% Unsaturated conductivity surface [mm/h]
        %         Lk = CF*Ks_b.*(((O-Ohy)./(Osat-Ohy)).^(1/2)).*...
        %             (1-(1-((O-Ohy)./(Osat-Ohy)).^(1./mvg)).^mvg).^2;
    else
        Lk = CF*Ks_b; %  [mm/h]
    end
    Lk(isnan(Lk)) = 0;
    
    leakage(i) = sum(Lk(:));
    % if leakage(i)>1000
    %     fprintf(':check:\n');
    % end
    Qi = kdtS.*To.*SLO/aT; % [mm/h] %%% Lateral Outflow
    
    V = V - Qi*dth - Lk*dth - ET*dth; % Volume second update [mm]
    
    Rd = V-Vmax; Rd(Rd<0)=0; % Dunne Runoff [mm]
    V = V-Rd; EX = -V.*(V<0); V(V<0)=0; % Volume Third update
    Rd = Rd/dth; % Dunne Runoff [mm/h]
    
    
    ex(i) = sum(EX(:));
    
    % Acquifer is for whole catchment.
    Vf(i) = Vf(i-1) + sum(sum(Lk*dth)) -sum(sum(EX)); %%% [mm] Total
%     if Vf(i)<0
%         disp('Time Step too Long')
%         break
%     end
    
    Qsur = Rd + Rh; %%% Surface flow [mm/h]
    
    % % Hide Permanantly @ yuting
    % CK1(i)= sum(sum(VI))/1000*(cellsize^2)  -  sum(sum(V))/1000*(cellsize^2) + sum(sum(Pr*dth))/1000*(cellsize^2) ...
    %     - sum(sum(ET*dth))/1000*(cellsize^2) + sum(sum(EX))/1000*(cellsize^2) ...
    %     - sum(sum(Qsur*dth))/1000*(cellsize^2) -  sum(sum(Qi*dth))/1000*(cellsize^2)...
    %     - sum(sum(Lk*dth))/1000*(cellsize^2);  %%[m^3]
    % CK1(i)=1000*CK1(i)/Area; %%[mm]
    % % Hide Permanantly @ yuting
    
    % Subsurface Routing 
    Qseep = Qi*dth.*(SN); %%[mm]
    Qi = Qi.*(1-SN); %% [mm/h]
    
    Qsurl = reshape(Qi*dth,DTM_nml,1);
    QiM = (XX)*Qsurl;
    QiM  = reshape(QiM , [mm nn]);
    
    QoutS(i) = (sum(sum(Qi*dth)) - sum(sum(QiM))); %%%  %%[mm]
    
    for p_i = 1:np
        QpointS(i,p_i) = Qi(Yout(p_i),Xout(p_i))*dth;%%%  %%[mm]
    end
    
    Qi_fall = 0*MASK;
    tH_store = 0*MASK;
    % OVERLAND FLOW
    Qsur = Qsur*dth; %%[mm]
    
    tc = 0;
    
    if any(Qsur(:))
        while tc<=dt
            
            tc = tc + dti;
            if tc > dt
                dti = dt-tc;
            end
            if dti > 0
                % SURFACE VELOCITY
                YH = (Qsur*dth/1000);% .*cellsize./WC ; %% [m]
                % Ponding - Microroughness
                
                YH = YH - MRO ; YH(YH<0)=0;
                
                tH = SR_m./(YH.^(2/3)); %%% [s]
                
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
                    Qpoint(i,p_i) = Qpoint(i,p_i) + kdt(Yout(p_i),Xout(p_i))*Qsur(Yout(p_i),Xout(p_i)); %%%%% -- [mm]
                end
                
                Qsur = QsurR; %[mm]
                Qsur = sparse(Qsur);
            end
            t1 = nanmin(nanmin(tH(SN==1)));
            dti = nanmin(nanmax(t1/5,2),60);
        end
    else
        YH = Qsur ; % [m]
        tH = Inf*Qsur; % [s]
        Qout(i) = 0; % -- [mm]
        Qpoint(i) = 0; %-- [mm]
        
        QsurR=0; %[mm]
        Qi_fall=0; %[mm]
    end
    
    
    Qch = QchR*dth + Qseep + Qi_fall; % [mm]
    % Qch = sparse(Qch);
    t_store(:) = 0;
    
    % Channel flow
    tc = 0;
    if any(Qch(:))
        
        while tc<=dt
            
            tc = tc+dti2;
            if tc>dt
                dti2 = dt-tc;
            end
            if dti2>0
                % SURFACE VELOCITY
                Y(SNl) = (Qch(SNl)).*CHR_m2(SNl) ; %% [m]
                % t(SNl)= (cellsize.*NMAN_C(SNl))./(Y(SNl).^(2/3).*SLO(SNl).^0.5); %%% [s]
                t(SNl) = CHR_m(SNl)./Y(SNl).^(2/3);
                t(isnan(t)) = 0;
                t_store = t_store + t*dti2/dt;
                kdt = dti2./t; %[]
                kdt(kdt>1) = 1;
                kdt(isnan(kdt)) = 0;
                
                Qchl = reshape(kdt.*Qch,DTM_nml,1);
                QchM = (XX)*(Qchl);
                QchM = reshape(QchM, [mm nn]);
                QchR = QchM + (Qch - kdt.*Qch) ; % [mm]
                
                QoutC(i) = QoutC(i) + (sum(sum(Qch))- sum(sum(QchR))); % -- [mm]
                for p_i = 1:np
                    QpointC(i,p_i) = QpointC(i,p_i) + kdt(Yout(p_i),Xout(p_i))*Qch(Yout(p_i),Xout(p_i)); %%%%% -- [mm]
                end
                Qch = QchR; %[mm]
                
                
            end
            t2 = nanmin(nanmin(t(SN == 1)));
            dti2 = nanmin(nanmax(t2/5,2),60);
        end
        
    else
        Y = Qch ; % [m]
        t = Qch*Inf; % [s]
        QoutC(i) = 0; % -- [mm]
        QpointC(i) = 0; % -- [mm]
        QchR = 0; %[mm]
    end
    
    % Show animation
    if tag_Animation
        spy(SN);
        hold on;
        imagesc(reshape(QpointC(i,:),56,68),[0,5]);
        colorbar;
        pause(0.001);
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
    O = V./Zs + Ohy; O(isnan(O))=0; % Soil Moisture update []
    QsurR = QsurR + Rd2/dth;  % [mm/h]
    
    if sum(isnan(QsurR(:)))>0
        disp(1);
    end
    
    % Acquifer Routing --- Linear
    % Qsub = Vf(i)/kF; % [mm/h] % simplified linear storage-discharge relationship
    
    % Acquifer Routing --- Non-Linear.
    Qsub = Vf(i).^(bgw)/kF;
    
    % Acquifer Routing --- Non-Linear with seasonality???
%     if Vf(i) > 5000%1000
%         Qsub = Vf(i).^(bgw)/kF;
%     else
%         kF1 = sum(SN(:))*[10.^1.15];
%         Qsub = Vf(i).^(1.15)/kF1;
%     end

    
    
    Vf(i) = Vf(i)- Qsub*dth; %[mm]
    QchR = QchR + SN*Qsub/sum(sum(SN));  % [mm/h]
    
    gwDisPart(i) = Qsub; % By Yuting
    
    % % Hide Permantly @ yuting
    % CK2(i)= sum(sum(VI))/1000*(cellsize^2)  -  sum(sum(V))/1000*(cellsize^2) + sum(sum(Pr*dth))/1000*(cellsize^2) ...
    %     - sum(sum(ET*dth))/1000*(cellsize^2) + sum(sum(EX))/1000*(cellsize^2)...
    %     -  sum(sum(QsurR*dth))/1000*(cellsize^2)  -  sum(sum(QchR*dth))/1000*(cellsize^2) ...
    %     - sum(sum(Lk*dth))/1000*(cellsize^2)+ Qsub/1000*(cellsize^2)*dth ...
    %     - QoutS(i)/1000*(cellsize^2)  -QoutC(i)/1000*(cellsize^2) -Qout(i)/1000*(cellsize^2) ;% - QpointS(i)/1000*(cellsize^2) -Qpoint(i)/1000*(cellsize^2);  %%[m^3]
    % CK2(i)=1000*CK2(i)/Area; %[mm]
    % % Hide Permantly @ yuting
    
    % Oave(i)= mean(reshape(O(O~=0),1,numel(O(O~=0))));
    % Ostd(i)= std(reshape(O(O~=0),1,numel(O(O~=0))));
    
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
    
    OUTPUT.SatSoil(i) = nanmean(reshape(catTag.*(O-Ohy)./(Osat-Ohy),[],1)>0.95);

    if pl == 1 & mod(i,100)==0
                
        % Plot soil moisture
        set(figure(1),'position',[50 -50 400 600])
        subplot(311)
        imagesc((O-Ohy)./(Osat-Ohy))
        caxis([0.2 1])
        axis off
        title(['Effective Saturation T = ' num2str(i) ' h'])
        % colorbar
        axis equal
        axis tight
        
        subplot(312)
        % plot((2:i),Qpoint(2:i)/1000*(cellsize^2)/dt,':m','Linewidth',1.5)
        % hold on ; grid on;
        % plot((2:i),QpointS(2:i)/1000*(cellsize^2)/dt,':g','Linewidth',1.5)
        % grid on ; hold on ;
        plot((2:i),QpointC(2:i)/1000*(cellsize^2)/dt,'r','Linewidth',1.5)
        grid on;
        ylabel(['Q [m^3/s]'])
        %ylim([0.001 1])
        
%         subplot(313)
%         plot((2:i),Vf(2:i),':m','Linewidth',1.5)
%         title(['Acquifer Storage']);
%         grid on ;
        
        pause(0.001)
    end
    
end

OUTPUT.gwDisPart = gwDisPart;
OUTPUT.leakage = leakage;
OUTPUT.Vf = Vf;
OUTPUT.ETr = ETr;
OUTPUT.QpointC = QpointC;
OUTPUT.Qpoint = Qpoint;
OUTPUT.QpointS = QpointS;
% OUTPUT.SatSoil has been saved;
% OUTPUT.O has been saved;
disp('FINE');
catch
1;    
end
end


