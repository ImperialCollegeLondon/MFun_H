clear all

load('PLM_TOPO_DATA')
sc = 0.25;

X_or = 0:15:(size(DEM_pl,1)-1)*15;
Y_or = 0:15:(size(DEM_pl,2)-1)*15;
X_or = X_or(end:-1:1);
[xx_or,yy_or] = meshgrid(Y_or,X_or);
F = scatteredInterpolant(xx_or(:),yy_or(:),DEM_pl(:));


% DEM_pl_c = imresize(DEM_pl,sc);
cellsize = 15*1/sc;
X = 0:cellsize:(size(DEM_pl,1)*sc-1)*cellsize;
X = X(end:-1:1);
Y = 0:cellsize:(size(DEM_pl,2)*sc-1)*cellsize;
[xx,yy] = meshgrid(Y,X);
DEM_pl_c = F(xx,yy);


%%




%%
DEM = GRIDobj(Y(:),X(:),DEM_pl_c);
FD  = FLOWobj(DEM,'preprocess','fill');
S   = STREAMobj(FD,'minarea',200);
DEM_f = fillsinks(DEM);
PIT = abs(DEM.Z-DEM_f.Z)>0;
[DTMF,FDR]= Flat_Area_Treatment(DEM.Z,PIT,0);
% %%
% IND = abs(DEM_f.Z-DEM.Z)>0;
% a = 1;
% DEM_f_t = DEM_f;
% while a
%     tic
%     DEM_f_t.Z(IND)=DEM_f_t.Z(IND)+1e-3*rand(sum(IND(:)),1);
%     DEM_f_t2 = fillsinks(DEM_f_t);
%     [R,S] = dem_flow(DEM_f_t.Z,cellsize,cellsize);
%     
%     if nansum( abs(DEM_f_t2.Z(:)-DEM_f_t.Z(:))>0 ) == 0 && sum(isnan(R(~isnan(DEM_f_t2.Z))))==0
%         a = 0;
%     end
%     DEM_f_t = DEM_f_t2;
%     toc
% end
%%
DEM_pl_c = DTMF;
DEM.Z = DEM_pl_c;
FD  = FLOWobj(DEM,'preprocess','fill');
% Use D-8
FA = flowacc(FD);
SN = double(log(FA.Z)>3.8);

DEM_pl_c(65,22) = 500;
[R,S] = dem_flow(DEM_pl_c,cellsize,cellsize);
[R,T] = D8_LTD(DEM_pl_c,cellsize,cellsize,0,'Transverse');
DEM_pl = DEM_pl_c;
x_pl = X;
y_pl = Y;
%%
save('PLM_TOPO_DATA_coarse','cellsize','DEM_pl','R','S','SN','T','x_pl','y_pl')