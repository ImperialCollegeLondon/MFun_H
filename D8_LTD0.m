function [R,T]=D8_LTD0(DTM,d1,d2,OPT)
%%%%%%%%%%%%%%%%
%%% Ref. Orlandini et al., 2003
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LTD with lambda = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n] = size(DTM);
[R,S] = dem_flow(DTM,d1,d2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dem_flow Downslope flow direction for a DEM
%   Reference: Tarboton, "A new method for the determination of flow
%   directions and upslope areas in grid digital elevation models," Water
%   Resources Research, vol. 33, no. 2, pages 309-319, February 1997.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p1=zeros(m,n);
p2=zeros(m,n);
sig=zeros(m,n);
del1=zeros(m,n);
del2=zeros(m,n);
alp1=zeros(m,n);
alp2=zeros(m,n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:(m)
    for j=1:(n)
        %%%%% Calculation angular deviation
        [p1(i,j),p2(i,j),sig(i,j),alp1(i,j),alp2(i,j)]=assign_AD(R(i,j));
        %%% Calculation transversal deviation
        del1(i,j)= d1*sin(alp1(i,j));
        del2(i,j)=sqrt(d1^2+d2^2)*sin(alp2(i,j));
        %%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Principal Directions
%%% p --> [ 1 2 3 4 5 6 7 8 9]
PriDir= [3*pi/4, pi, 5*pi/4,  pi/2, -1,  3*pi/2,  pi/4, 0, 7*pi/4];
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:(m)
    for j=1:(n)
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        if sum(ismember(R(i,j),PriDir))==0 && not(isnan(R(i,j)))
            %%%%%%%%%%%%%%%%%%%%%%%%
            switch OPT
                case 'Transverse'
                    if abs(del1(i,j)) <=  abs(del2(i,j))
                        R(i,j)=PriDir(p1(i,j));
                    else
                        R(i,j)=PriDir(p2(i,j));
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%
                case 'Angular'
                    if abs(alp1(i,j)) <=  abs(alp2(i,j))
                        R(i,j)=PriDir(p1(i,j));
                    else
                        R(i,j)=PriDir(p2(i,j));
                    end
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = flow_matrix(DTM,R,d1,d2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p1,p2,sig,alp1,alp2]=assign_AD(R)
if R >=0 && R < pi/4
    p1=8;
    p2=7;
    sig=-1;
    alp1=R;
    alp2=pi/4-R;
elseif  R >=pi/4 && R < pi/2
    p1=4;
    p2=7;
    sig=1;
    alp2=R-pi/4;
    alp1=pi/2-R;
elseif   R >=pi/2 && R < 3*pi/4
    p1=4;
    p2=1;
    sig=-1;
    alp1=R-pi/2;
    alp2=3*pi/4-R;
elseif   R >=3*pi/4 && R < pi
    p1=2;
    p2=1;
    sig=1;
    alp2=R-3*pi/4;
    alp1=pi-R;
elseif   R >=pi && R < 5*pi/4
    p1=2;
    p2=3;
    sig=-1;
    alp1=R-pi;
    alp2=5*pi/4-R;
elseif  R >=5*pi/4 && R < 3*pi/2
    p1=6;
    p2=3;
    sig=1;
    alp2=R-5*pi/4;
    alp1=3*pi/2-R;
elseif   R >=3*pi/2 && R < 7*pi/4
    p1=6;
    p2=9;
    sig=-1;
    alp1=R-3*pi/2;
    alp2=7*pi/4-R;
elseif   R >=7*pi/4 && R < 2*pi
    p1=8;
    p2=9;
    sig=1;
    alp2=R-7*pi/4;
    alp1=2*pi-R;
elseif isnan(R)
    p1=NaN;
    p2=NaN;
    sig=NaN;
    alp2=NaN;
    alp1=NaN;
end
end


