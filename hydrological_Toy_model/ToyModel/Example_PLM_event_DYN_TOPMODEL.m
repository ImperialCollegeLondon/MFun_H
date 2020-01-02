close all
clear all

% Load the data and set up model parameters

load PLM_TOPO_DATA2.mat;
load PLM_data_h;
SIM_n = 5*24;

DTM = DEM_pl;
%%% Precipitation
Prt = zeros(SIM_n,1);
Prt(5:10) = 10;
Prt(48:53) = 10;

%%% Evapotranspiration
ETPt = PET_h(1:SIM_n);

[m,n] = size(DTM);
MASK=ones(size(DTM)); 
MASK(isnan(DTM))=0;
SN(isnan(SN))=0;

Zs = zeros(m,n);
Zs(:,:) = 400.*MASK; %%  Soil Active depth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Osat =  (zeros(m,n) + 0.35).*MASK;
Ohy =   (zeros(m,n) + 0.1).*MASK;
Oel =   (zeros(m,n) + 0.11).*MASK;
%%%%%%%%%%%%%%%%%%
nvg= 2.8;
Ks=  (zeros(m,n) + 20).*MASK; %%% [mm/h]
mvg = (zeros(m,n) + (1 - 1/nvg)).*MASK;
Kbot=  (zeros(m,n) + 0.01).*MASK; %%% [mm/h]
%%%%%%%%%%%%%%%%%%%%%%%%%%
kF = 10000;%   %%% [h]  acquifer constant
mF= 100000; %320; %%% [mm]
%%%%%
OPT_UNSAT = 1; %%% Option unsaturated 1 -- saturated bottom 0
CF=0; %% Acquifer Yes - No
aR=100; %%[-] anisotropy ratio

dt =3600;%%%[s]
dti= 10; %%[s] Internal Time step for Surface Channel-flow Routing
dti2= 2; %%[s] Internal Time step for Surface Overland-flow Routing
T_map_s = 10*24;
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

Xout=205; Yout=277;

%%%%%%% Intial Oi
Oi=0.1.*MASK; 
O=Oi;

%% SET UP SIMULATION

TOPO_DAT.DTM = DTM;
TOPO_DAT.T = T;
TOPO_DAT.SN = SN;
TOPO_DAT.WC = WC;
TOPO_DAT.MASK = MASK;
TOPO_DAT.Zs = Zs; 
TOPO_DAT.cellsize = cellsize;

HYDR_DATA.Osat = Osat;
HYDR_DATA.Ohy = Ohy;
HYDR_DATA.Oel = Oel;
HYDR_DATA.Ks = Ks;
HYDR_DATA.aR = aR; %%[-] anisotropy ratio
HYDR_DATA.mvg = mvg;
HYDR_DATA.Kbot = Kbot;
HYDR_DATA.kF = kF;%   %%% [h]  acquifer constant
HYDR_DATA.mF = mF; %320; %%% [mm]
HYDR_DATA.NMAN_C = NMAN_C; 
HYDR_DATA.NMAN_H = NMAN_H; 
HYDR_DATA.MRO = MRO; %0.002; %%[m]

METEO_DAT.Prt = Prt;
METEO_DAT.ETPt = ETPt;
METEO_DAT.D_h = D_h;

SIM_PARAM.dt = dt;
SIM_PARAM.dti = dti;
SIM_PARAM.dti2 = dti2;
SIM_PARAM.OPT_UNSAT = OPT_UNSAT;
SIM_PARAM.CF = CF; %% Acquifer Yes - No
SIM_PARAM.Oi = Oi;
SIM_PARAM.pl = 1;

OUT_PARAM.Xout = Xout;
OUT_PARAM.Yout = Yout;
OUT_PARAM.T_map_s = T_map_s;


%% RUN SIMULATION 
% [OUTPUT] = DYN_TOPMODEL(TOPO_DAT,HYDR_DATA,METEO_DAT,SIM_PARAM,OUT_PARAM);
 [OUTPUT] = DYN_TOPMODEL_ad_ts(TOPO_DAT,HYDR_DATA,METEO_DAT,SIM_PARAM,OUT_PARAM);

%% Plot output

figure(1)
ZZ = OUTPUT.QpointC(2:end)/1000*(cellsize^2)/dt;
% ZZ = Qpoint(2:end)/1000*(cellsize^2)/dt;
plot([Runoff(1:ceil(SIM_n/24))],'kx')
hold on
plot(mean(buffer(ZZ(:),24)),'r-')
xlabel('Time [d]'); ylabel('Discharge [m^3/s]');
