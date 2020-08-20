%% Distributed Data sets
% Tutorial week 6

% add the TopoToolbox
addpath(genpath('topotoolbox'));

%% Read Spatial Data

% DEM = GRIDobj('DEM.asc');
% [X,Y] = refmat2XY_p(DEM.refmat,DEM.size); 
% [xx,yy] = meshgrid(X,Y);
D = load('PLM_TOPO_DATA');
DEM = GRIDobj(D.x_pl(1,:),D.y_pl(:,1),D.DEM_pl);

xx = D.x_pl;
yy = D.y_pl;
X = D.x_pl(1,:);
Y = D.y_pl(:,1);
FD  = FLOWobj(DEM,'preprocess','fill');
S   = STREAMobj(FD,'minarea',200);


imageschs(DEM)

%% Exploratory plots

% Elevation plot
figure(1)
h = pcolor(X,Y,DEM.Z);
set(h, 'EdgeColor', 'none');
xlabel('X [m]')
ylabel('Y [m]')
title('Elevation')
colorbar
axis equal 
axis tight

% Gradient
Z_gr = gradient8(DEM,'deg');
figure(2)
h = pcolor(X,Y,Z_gr.Z);
set(h, 'EdgeColor', 'none');
title('Steepest gradient')
xlabel('X [m]')
ylabel('Y [m]')
colorbar
axis equal 
axis tight

% Aspect
Asp = aspect(DEM);
figure(3)
h = pcolor(X,Y,Asp.Z);
set(h, 'EdgeColor', 'none');
title('Aspect')
xlabel('X [m]')
ylabel('Y [m]')
colorbar
axis equal 
axis tight

% Curvature
Cur = curvature(DEM);
figure(31)
h = pcolor(X,Y,Cur.Z);
set(h, 'EdgeColor', 'none');
title('Curvature')
xlabel('X [m]')
ylabel('Y [m]')
caxis([-0.01 0.01])
colorbar
axis equal 
axis tight
%% Fill Pits
DEM = fillsinks(DEM);

%% Accumulation areas
FD  = FLOWobj(DEM,'preprocess','fill');
% Use D-8
FA = flowacc(FD);

figure(4)
subplot(121)
h = pcolor(X,Y,log(FA.Z));
set(h, 'EdgeColor', 'none');
title('D-8')
xlabel('X [m]')
ylabel('Y [m]')
colorbar
axis equal 
axis tight

% Use D-Inf
FD  = FLOWobj(DEM,'Dinf');
FA = flowacc(FD);

figure(4)
subplot(122)
h = pcolor(X,Y,log(FA.Z));
set(h, 'EdgeColor', 'none');
title('D-\infty')
xlabel('X [m]')
ylabel('Y [m]')
colorbar
axis equal 
axis tight

% Strahler order
figure(5)
FD  = FLOWobj(DEM,'preprocess','fill');
FA = flowacc(FD)*FD.cellsize^2;
W = FA.Z>45000;
[S] = streamorder(FD,W);
h = pcolor(X,Y,S.Z);
set(h, 'EdgeColor', 'none');
title('Strahler order')
xlabel('X [m]')
ylabel('Y [m]')
colorbar
axis equal 
axis tight

% Flow Distance
Flow_d = flowdistance(FD);
figure(6)
h = pcolor(X,Y,Flow_d.Z);
set(h, 'EdgeColor', 'none');
title('Flow Distance')
xlabel('X [m]')
ylabel('Y [m]')
colorbar
axis equal 
axis tight

% Catchment dileneation
figure(61)
X_cor = 282912;
Y_cor = 283637;
[L,outlet] = drainagebasins(FD,X_cor,Y_cor);
h = pcolor(X,Y,L.Z);
set(h, 'EdgeColor', 'none');
hold on
contour(X,Y,DEM.Z,20,'w')
plot(xx(outlet),yy(outlet),'ko','markerfacecolor','y')
caxis([0 max(L.Z(:))+1])
xlabel('X [m]')
ylabel('Y [m]')
axis equal 
axis tight
%% Hillshading effects

% Sky view factor
[HZ,Za] = Horizon_Angle(DEM.Z,DEM.cellsize);
[SvF,Ct] = Sky_View_Factor(DEM.Z,Z_gr.Z,Asp.Z,HZ,Za);

figure(7)
h = pcolor(X,Y,SvF);
set(h, 'EdgeColor', 'none');
title('Sky view factor')
xlabel('X [m]')
ylabel('Y [m]')
caxis([0.7 1])
colorbar
axis equal 
axis tight

% Shadow area 
Datam = [2017 06 22 20];
DeltaGMT = 0;
Lat = 52.5;
Lon = -3;
t_bef = 0.5;
t_aft = 0.5;
[h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day] = SetSunVariables(Datam,DeltaGMT,Lon,Lat,t_bef,t_aft);


[ShF] = Shadow_Effect(DEM.Z,h_S,zeta_S,HZ,Za);
figure(8)
h = pcolor(X,Y,ShF);
set(h, 'EdgeColor', 'none');
title('Shadow')
xlabel('X [m]')
ylabel('Y [m]')
caxis([0.0 1])
colorbar
axis equal 
axis tight

%% Clear Sky Radiation throughout a day

Z_gr = gradient8(DEM,'tan');
Slo_top = Z_gr.Z;
Aspect = Asp.Z;

DN = datenum([2018 6 21 0 0 0]):1/24:datenum([2018 6 21 23 0 0]);
Zbas = nanmean(DEM.Z(:));
Tdew = 20;
So1 = 1366*0.465;
So2 = 1366*0.52;
beta_A = 0.05;
alpha_A  = 1.3;
omega_A1  = 0.92;
omega_A2 = 0.84;
uo = 0.35;
un  = 0.0002;
rho_g = 0.15;
    
for i = 1:length(DN)
    
    Datam = datevec(DN(i));
    Datam = Datam(1:4);
    
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day] = SetSunVariables(Datam,DeltaGMT,Lon,Lat,t_bef,t_aft);
    [ShF] = Shadow_Effect(DEM.Z,h_S,zeta_S,HZ,Za);

    [Eb1,Eb2,Ed1,Ed2,Edp1,Edp2,rho_s1,rho_s2,Mb,Mg] ...
        = SetClearSkyRadiation(h_S,Zbas,Tdew,So1,So2,beta_A,alpha_A,omega_A1,omega_A2,uo,un,rho_g);
    
    cos_fst = cos(atan(Slo_top))*sin(h_S) + sin(atan(Slo_top)).*cos(h_S).*cos(zeta_S-Aspect*pi/180);

    if sin(h_S) > 0.10
        SAD1_S = Ed1*SvF + Ct.*rho_g.*((Ed1/sin(h_S)).*cos_fst + (1-SvF).*Ed1);
        SAD2_S = Ed2*SvF + Ct.*rho_g.*((Ed2/sin(h_S)).*cos_fst + (1-SvF).*Ed2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SAB1_S =(Eb1/sin(h_S)).*cos_fst.*ShF;
        SAB2_S =(Eb2/sin(h_S)).*cos_fst.*ShF;
    else
        SAD1_S = zeros(size(DEM.Z));
        SAD2_S = zeros(size(DEM.Z));
        SAB1_S = zeros(size(DEM.Z));
        SAB2_S = zeros(size(DEM.Z));
    end
    
    set(figure(10),'position',[10 10 1500 800])
    h=surf(X,Y,DEM.Z,SAD1_S+SAD2_S+SAB1_S+SAB2_S);
    set(h, 'EdgeColor', 'none');
    title(['I [W/m^2] H= ' num2str(Datam(4))])
    xlabel('X [m]')
    ylabel('Y [m]')
    caxis([0 1300])
    colorbar
    view(30,30)
    pause
    
end
