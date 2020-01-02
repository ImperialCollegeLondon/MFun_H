function [R2,T]=D8_LTD(DTM,d1,d2,lambda,OPT)
%%%%%%%%%%%%%%%%
%%% Ref. Orlandini et al., 2003; Orlandini and Moretti 2009 
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LTD with lambda = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n]=size(DTM);
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
deccum=zeros(m,n);
w1=zeros(m,n);
w2=zeros(m,n);
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
ich = [ -1 0 1 -1 0 1 -1 0 1];
jch = [ -1 -1 -1 0 0 0 +1 +1 +1];
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
switch OPT
    case 'Transverse'
        %lambda = 1;
        [qta,Ind]=sort(reshape(-DTM,1,m*n));
        Ind=Ind(not(isnan(qta)));
        R2=R*NaN;
        FACC=DTM*0;
        T=speye(m*n,m*n);
        bau = waitbar(0,'Internal Loop');
        for k=1:1:(length(Ind)-1)
            waitbar(k/length(Ind),bau);
            %%%%%%%%%%%
            if not(isnan(DTM(Ind(k))))
                [i,j]= ind2sub([m,n],Ind(k));
                FACC(i,j)=FACC(i,j)+1;
                %%%%%%%%%%%%%%%%%%%%%%%%%
                if  FACC(i,j) == 1
                    %%%%%%%%%%%%%%%%%%%%%%%%
                    dec1 = sig(i,j)*del1(i,j);
                    dec2 = -sig(i,j)*del2(i,j);
                    %%%%%%%%%%%%%%%%%%%%%
                    if abs(dec1) <=  abs(dec2)
                        R2(i,j)=PriDir(p1(i,j));
                        deccum(i,j) = dec1;
                        FACC(i+ich(p1(i,j)),j+jch(p1(i,j)))= FACC(i+ich(p1(i,j)),j+jch(p1(i,j)))+FACC(i,j);
                        njk=sub2ind([m,n],i+ich(p1(i,j)),j+jch(p1(i,j)));
                        T(njk,Ind(k))=-1;
                    else
                        R2(i,j)=PriDir(p2(i,j));
                        deccum(i,j) = dec2;
                        FACC(i+ich(p2(i,j)),j+jch(p2(i,j)))= FACC(i+ich(p2(i,j)),j+jch(p2(i,j)))+FACC(i,j);
                        njk=sub2ind([m,n],i+ich(p2(i,j)),j+jch(p2(i,j)));
                        T(njk,Ind(k))=-1;
                    end
                else
                    Ic = find(T(Ind(k),:)~=0);
                    Ic = setdiff(Ic,Ind(k));
                    if isempty(Ic);
                        disp('Error no upstream contributions')
                        return ;
                    end
                    [iv,jv]= ind2sub([m,n],Ic);
                    %%%%%%%%%%%%%%%%%%%%%%%%
                    %dec1 = sig(i,j)*del1(i,j) + sum(lambda*FACC(iv,jv).*deccum(iv,jv)/sum(FACC(iv,jv)));
                    %dec2 = -sig(i,j)*del2(i,j)+ sum(lambda*FACC(iv,jv).*deccum(iv,jv)/sum(FACC(iv,jv)));
                    dec1 = sig(i,j)*del1(i,j) + lambda*sum(FACC(iv,jv).*deccum(iv,jv))/sum(FACC(iv,jv));
                    dec2 = -sig(i,j)*del2(i,j)+ lambda*sum(FACC(iv,jv).*deccum(iv,jv))/sum(FACC(iv,jv));
                    %%%%%%%%%%%%%%%%%%%%%
                    if abs(dec1) <=  abs(dec2)
                        R2(i,j)=PriDir(p1(i,j));
                        deccum(i,j) = dec1;
                        FACC(i+ich(p1(i,j)),j+jch(p1(i,j)))= FACC(i+ich(p1(i,j)),j+jch(p1(i,j)))+FACC(i,j);
                        njk=sub2ind([m,n],i+ich(p1(i,j)),j+jch(p1(i,j)));
                        T(njk,Ind(k))=-1;
                    else
                        R2(i,j)=PriDir(p2(i,j));
                        deccum(i,j) = dec2;
                        FACC(i+ich(p2(i,j)),j+jch(p2(i,j)))= FACC(i+ich(p2(i,j)),j+jch(p2(i,j)))+FACC(i,j);
                        njk=sub2ind([m,n],i+ich(p2(i,j)),j+jch(p2(i,j)));
                        T(njk,Ind(k))=-1;
                    end
                end
            end
        end
        close(bau);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Angular'
        %lambda = 1;
        [qta,Ind]=sort(reshape(-DTM,1,m*n));
        Ind=Ind(not(isnan(qta)));
        R2=R*NaN;
        FACC=DTM*0;
        T=speye(m*n,m*n);
        bau = waitbar(0,'Internal Loop');
        for k=1:1:(length(Ind)-1)
            waitbar(k/length(Ind),bau);
            %%%%%%%%%%%
            if not(isnan(DTM(Ind(k))))
                [i,j]= ind2sub([m,n],Ind(k));
                FACC(i,j)=FACC(i,j)+1;
                %%%%%%%%%%%%%%%%%%%%%%%%%
                if  FACC(i,j) == 1
                    %%%%%%%%%%%%%%%%%%%%%%%%
                    dec1 = sig(i,j)*alp1(i,j);
                    dec2 = -sig(i,j)*alp2(i,j);
                    %%%%%%%%%%%%%%%%%%%%%
                    if abs(dec1) <=  abs(dec2)
                        R2(i,j)=PriDir(p1(i,j));
                        deccum(i,j) = dec1;
                        FACC(i+ich(p1(i,j)),j+jch(p1(i,j)))= FACC(i+ich(p1(i,j)),j+jch(p1(i,j)))+FACC(i,j);
                        njk=sub2ind([m,n],i+ich(p1(i,j)),j+jch(p1(i,j)));
                        T(njk,Ind(k))=-1;
                    else
                        R2(i,j)=PriDir(p2(i,j));
                        deccum(i,j) = dec2;
                        FACC(i+ich(p2(i,j)),j+jch(p2(i,j)))= FACC(i+ich(p2(i,j)),j+jch(p2(i,j)))+FACC(i,j);
                        njk=sub2ind([m,n],i+ich(p2(i,j)),j+jch(p2(i,j)));
                        T(njk,Ind(k))=-1;
                    end
                else
                    Ic = find(T(Ind(k),:)~=0);
                    Ic = setdiff(Ic,Ind(k));
                    if isempty(Ic);
                        disp('Error no upstream contributions')
                        return ;
                    end
                    [iv,jv]= ind2sub([m,n],Ic);
                    %%%%%%%%%%%%%%%%%%%%%%%%
                    dec1 = sig(i,j)*alp1(i,j) + lambda*sum(FACC(iv,jv).*deccum(iv,jv))/sum(FACC(iv,jv));
                    dec2 = -sig(i,j)*alp2(i,j)+ lambda*sum(FACC(iv,jv).*deccum(iv,jv))/sum(FACC(iv,jv));
                    %%%%%%%%%%%%%%%%%%%%%
                    if abs(dec1) <=  abs(dec2)
                        R2(i,j)=PriDir(p1(i,j));
                        deccum(i,j) = dec1;
                        FACC(i+ich(p1(i,j)),j+jch(p1(i,j)))= FACC(i+ich(p1(i,j)),j+jch(p1(i,j)))+FACC(i,j);
                        njk=sub2ind([m,n],i+ich(p1(i,j)),j+jch(p1(i,j)));
                        T(njk,Ind(k))=-1;
                    else
                        R2(i,j)=PriDir(p2(i,j));
                        deccum(i,j) = dec2;
                        FACC(i+ich(p2(i,j)),j+jch(p2(i,j)))= FACC(i+ich(p2(i,j)),j+jch(p2(i,j)))+FACC(i,j);
                        njk=sub2ind([m,n],i+ich(p2(i,j)),j+jch(p2(i,j)));
                        T(njk,Ind(k))=-1;
                    end
                end
            end
        end
        close(bau);
end

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


