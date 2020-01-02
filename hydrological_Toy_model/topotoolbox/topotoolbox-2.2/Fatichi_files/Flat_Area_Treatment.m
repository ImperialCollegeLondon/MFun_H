function[DTMF,FDR]= Flat_Area_Treatment(DTM,PIT,opt_out_bou)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% It works for DTM with elevation z > 0 m
%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n]=size(DTM);
DTMf= [ NaN*ones(1,n+2); [NaN*ones(m,1), DTM, NaN*ones(m,1)] ; NaN*ones(1,n+2) ];
%%%%%%% Option outlet in the boundary %%%%%%%%%%%%%%%%%%%%%
if opt_out_bou == 1
    TEM = DTMf>0;
    BOU= circshift(TEM,[0,1]) +  circshift(TEM,[-1,1])+  circshift(TEM,[-1,0])+ ...
        circshift(TEM,[-1,-1])+  circshift(TEM,[0,-1])+ ...
        circshift(TEM,[1,-1])+ circshift(TEM,[1,0]) + circshift(TEM,[1,1]);
    BOU(BOU>1)=1; BOU = BOU-TEM;
    BOU(BOU>1)=1; BOU(BOU<0)=0;
    DTMf(BOU==1)= nanmin(nanmin(DTMf))-1;
    clear BOU TEM
end
PIT= [ NaN*ones(1,n+2); [NaN*ones(m,1), PIT, NaN*ones(m,1)] ; NaN*ones(1,n+2) ];
CK = PIT;
FDR =CK*NaN;
DTMF = DTMf;
fdr_p = [ 128 64 32 ; 1 -9999 16 ; 2 4 8 ];
for i=2:(m+1)
    for j=2:(n+1)
        if not(isnan(DTMf(i,j))) &&  (CK(i,j) == 1)
            stoop = 0;
            iP = [i];
            jP=  [j];
            CK(i,j)=0;
            Pex= CK(i-1:i+1,j-1:j+1);
            f_ind=find(Pex==1);
            liP_old = length(iP); i2 = i; j2= j;
            while stoop == 0
                [iT,jT] = ind2sub([3,3],f_ind);
                for r=1:length(iT)
                    iP=[iP i2+iT(r)-2];
                    jP=[jP j2+jT(r)-2];
                    CK(i2+iT(r)-2,j2+jT(r)-2)= 0;
                end
                iT=[]; jT=[];
                if length(iP)> liP_old
                    i2= iP(liP_old + 1);
                    j2= jP(liP_old + 1);
                    Pex= CK(i2-1:i2+1,j2-1:j2+1);
                    f_ind=find(Pex==1);
                    liP_old = liP_old + 1;
                else
                    stoop = 1 ;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            K = sub2ind(size(PIT),iP,jP);
            TEM = zeros(m+2,n+2); TEM(K)=1;
            BOU= circshift(TEM,[0,1]) +  circshift(TEM,[-1,1])+  circshift(TEM,[-1,0])+ ...
                circshift(TEM,[-1,-1])+  circshift(TEM,[0,-1])+ ...
                circshift(TEM,[1,-1])+ circshift(TEM,[1,0]) + circshift(TEM,[1,1]);
            BOU(BOU>1)=1; BOU = BOU-TEM; 
            BOU(BOU>1)=1; BOU(BOU<0)=0; 
            Z = DTMf.*BOU;  Z(Z==0)=NaN;
            [Zmin, p1] = nanmin(Z); [Zmin, p2] = nanmin(Zmin);
            KP=sub2ind(size(PIT),p1(p2),p2);
            Z(Z==Zmin)=Inf;
            Zmin2 = nanmin(nanmin(Z));
            DZ = Zmin2-Zmin; DZ=DZ/(length(K)+2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            k=1; 
            while not(isempty(K))
                for r=1:length(K)
                    switch K(r)
                        case KP(k)-1
                            FDR(K(r))=fdr_p(4);
                            Zmin=Zmin+DZ;
                            DTMF(K(r))=Zmin; KP=[KP K(r)]; 
                        case KP(k)+1
                            FDR(K(r))=fdr_p(6);
                            Zmin=Zmin+DZ;
                            DTMF(K(r))=Zmin; KP=[KP K(r)]; 
                        case KP(k) - (m+2)
                            FDR(K(r))=fdr_p(2);
                            Zmin=Zmin+DZ;
                            DTMF(K(r))=Zmin; KP=[KP K(r)]; 
                        case KP(k) - (m+2)-1
                            FDR(K(r))=fdr_p(1);
                            Zmin=Zmin+DZ; 
                            DTMF(K(r))=Zmin; KP=[KP K(r)]; 
                        case KP(k) - (m+2)+1
                            FDR(K(r))=fdr_p(3);
                            Zmin=Zmin+DZ;
                            DTMF(K(r))=Zmin;KP=[KP K(r)]; 
                        case KP(k) + (m+2)
                            FDR(K(r))=fdr_p(8);
                            Zmin=Zmin+DZ;
                            DTMF(K(r))=Zmin;KP=[KP K(r)]; 
                        case KP(k) + (m+2)-1
                            FDR(K(r))=fdr_p(7);
                            Zmin=Zmin+DZ;
                            DTMF(K(r))=Zmin;KP=[KP K(r)]; 
                        case KP(k) + (m+2)+1
                            FDR(K(r))=fdr_p(9);
                            Zmin=Zmin+DZ;
                            DTMF(K(r))=Zmin;KP=[KP K(r)]; 
                    end
                end 
                for r=1:length(KP); 
                    K=setdiff(K,KP(r));
                end
                k=k+1; 
            end
        end
        clear K TEM BOU KP Z Zmin Zmin2 DZ  iP jP Pex liP_old f_ind k r 
    end
end
FDR=FDR(2:end-1,2:end-1);
DTMF=DTMF(2:end-1,2:end-1);
end
