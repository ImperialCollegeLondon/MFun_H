% Example_PLM_DYN_TOPMODEL.m provides a simplified model to compute the
% discharge using the grid2grid hydrological model.

close all
clear all

% Load the data and set up model parameters

load PLM_TOPO_DATA_coarse.mat; 
load PLM_data_h;
SIM_s = 6500;% starting index of used data;
SIM_n = 1000;% simulation length;

DTM = DEM_pl;
%%% Precipitation
Prt = Prs(SIM_s:SIM_s+SIM_n); %%% [mm/h]

%%% Evapotranspiration
ETPt = PET_h(SIM_s:SIM_s+SIM_n);

[m,n] = size(DTM);
MASK = ones(size(DTM)); 
MASK(isnan(DTM)) = 0;
SN(isnan(SN)) = 0;%Stream Network

Zs = zeros(m,n);
Zs(:,:) = 400.*MASK;%%  Soil Active depth; [mm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%% <van Genuchten related>
Osat = (zeros(m,n) + 0.35).*MASK;% Saturated Water Content of Soil.
Ohy =  (zeros(m,n) + 0.1).*MASK;% Residual Water Content (WP Water Content)
Oel =  (zeros(m,n) + 0.11).*MASK;%<(zeros(m,n) + 0.11).*MASK;>% FC Water Content
Psi_ae = (zeros(m,n) + 0.11).*MASK;% Suction pressure
% @ Yuting
% Parameters for more soil textural type can be found in 'Table 2'.
% Availabel from: https://ac.els-cdn.com/S0168192314000483/1-s2.0-S0168192314000483-main.pdf?
%    _tid=cbbb639b-ed8d-4616-83de-8177a582cb65&acdnat=1548284892_b51a1a93870c46b902286564ba17d5d9
% Reference: Modeling plant transpiration under limited soil water: Comparison of different plant and 
%    soil hydraulic parameterizations and preliminary implications for their use in land surface models

%%%%%%%%%%%%%%%%%%%%%%%%%%%% <Soil Conductivity related>
nvg= 2.8;
Ks=  (zeros(m,n) + 20).*MASK; %%% [mm/h] % Saturated Hydraulic Conductivity;
mvg = (zeros(m,n) + (1 - 1/nvg)).*MASK;
Kbot=  (zeros(m,n) + 0.01).*MASK; %%% [mm/h]

%%%%%%%%%%%%%%%%%%%%%%%%%%%% <Acquifer Conductivity related>
kF = 10000; %%% [h]  acquifer constant
mF= 100000; %320; %%% [mm] ????
%%%%%%%%%%%%%%%%%%%%%%%%%%%% <Routing related>
OPT_UNSAT = 1; %%% Option unsaturated 1 -- saturated bottom 0
CF = 0; %% Acquifer Yes - No
aR=100; %%[-] anisotropy ratio

dt = 3600;%%[s]
dti = 10; %%[s] Internal Time step for Surface Channel-flow Routing
dti2 = 2; %%[s] Internal Time step for Surface Overland-flow Routing
T_map_s = 10*24;% 10 days?
GA = 1;
%%%%%%%%%%%%%%%%%%%%%%%%

DEM = GRIDobj(y_pl(:),x_pl(:),DEM_pl);
FD  = FLOWobj(DEM);
FA = flowacc(FD);
WC = zeros(m,n); %%% [m]  width channel
WC(SN==1) = 0.0015*sqrt(FA.Z(SN==1)*cellsize^2); %% [m]

%%%%%%%%%%%%%%%%%%%%%%%%
WC = WC.*MASK;
NMAN_C = SN.*0.03;
NMAN_H = 0.10;
NMAN_C = NMAN_C.*MASK;%%[s/(m^1/3)] manning coefficient
NMAN_H = NMAN_H.*MASK;%%???
MRO = 0.0; %0.002; %%[m] Manning Roughness                            

Xout = 51; Yout = 68;

%%%%%%% Intial Oi
Oi = 0.26.*MASK;% initial soil moisture;                    
O = Oi;

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
HYDR_DATA.Psi_ae = 0.5;
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
SIM_PARAM.pl = 0;
SIM_PARAM.GA = GA;

OUT_PARAM.Xout = Xout;% Xout=51; Yout=68;
OUT_PARAM.Yout = Yout;
OUT_PARAM.T_map_s = T_map_s;


%% RUN SIMULATION 
% [OUTPUT] = DYN_TOPMODEL(TOPO_DAT,HYDR_DATA,METEO_DAT,SIM_PARAM,OUT_PARAM);
 [OUTPUT] = DYN_TOPMODEL_ad_ts(TOPO_DAT,HYDR_DATA,METEO_DAT,SIM_PARAM,OUT_PARAM);

%% Plot output
figure(1);
pointi = 1;
ZZ = OUTPUT.QpointC(2:end,pointi)/1000*(cellsize^2)/dt;
% ZZ = Qpoint(2:end)/1000*(cellsize^2)/dt;
plot([Runoff(ceil(SIM_s/24):ceil((SIM_s+SIM_n)/24))],'ko:')
hold on
plot(mean(buffer(ZZ(:),24)),'r-')
xlabel('Time [d]'); ylabel('Discharge [m^3/s]');
legend('PLM Observed','Simulated');
