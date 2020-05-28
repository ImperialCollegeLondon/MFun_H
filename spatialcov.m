%% This function returns the spatial covariance of the signal val with 
% locations X and Y
% The covariance is estimated at binned intervals with width ds out to a
% maximum distance smax
% The outputs of the function are the spatial covariance values and the
% covariance distance.
function [covariance,distance,r,disI,disJ,Ind_p1s,Ind_p2s]=spatialcov(X,Y,Data,ds,smax,r,disI,disJ,Ind_p1s,Ind_p2s)
%
% Note:
% This script is designed for QUICKLY computing spatial autocovariance from a
% huge amount of images.
% @ Yuting
%

if length(size(Data)) == 2
    Data = reshape(Data,[size(Data),1]);
end
if prod(size(Data,[1,2]))>12100 || prod(size(Data,[1,2,3]))>12100*20000
    error('Normal PCs Memory might be not enough to compute this BIG array in this fast-script.')
end

X = X(:);
Y = Y(:);
val = Data(:);
vest=std(val);
disp(['Mean: ',num2str(mean(val)),' Variance: ',num2str((vest))])
covariance = zeros(round(smax/ds)+2,size(Data,3));
ncov = zeros(round(smax/ds)+2,size(Data,3));
distance = [1:smax]';
if isempty(r)
    r = zeros(length(X),length(X));
    [disI,disJ] = deal(cell(1+smax,1));
    for p1_i = 1:length(X)
        for p2_i = 1:length(X)
            if p1_i <= p2_i % prevent double counting
                r(p1_i,p2_i) = sqrt((X(p1_i)-X(p2_i)).^2+(Y(p1_i)-Y(p2_i)).^2);
            end
        end
    end
    r = int16(round(r));
    [Ind_p1s,Ind_p2s] =  deal([]);
    for ir = 1:length(distance)
        [p1_s,p2_s] = find(r==distance(ir));
        Ind_p1s{ir,1} = int16(p1_s);
        Ind_p2s{ir,1} = int16(p2_s);
    end
end
Data = reshape(Data,[],size(Data,3));
for ir = 1:length(distance)
    tic
    p1_s = Ind_p1s{ir};
    p2_s = Ind_p2s{ir};
    thisWholeLen = length(p1_s);
    setStart = 1;
    setLen = 500000;% this number is set based on yt PC's total RAM.
    for setNo = 1:ceil(thisWholeLen/setLen)
        setInd = setStart:1:min(setStart+setLen-1,thisWholeLen);
        if ~isempty(setInd)
            covariance(ir,:) = covariance(ir,:)+...
                sum(Data(p1_s(setInd),:)&Data(p2_s(setInd),:),1);
            setStart = setInd(end)+1;
        end
    end
    ncov(ir,:) = ncov(ir,:)+thisWholeLen;
    toc
end

% for i = 1:size(D,1)
%     for j = 1:size(D,1)
%         if r(i,j)<smax && j~=i
%             ir = round((r(i,j)/ds))+1;
%             covariance(ir,:) = covariance(ir,:) + D(i,:) & D(j,:) ;
%             ncov(ir,:) = ncov(ir,:)+1;
%         end
%     end
% end

covariance = covariance./ncov;  
covariance((length(distance)+1):end,:) = [];

% figure
plot(distance,nanmean(covariance,2),'o')
xlabel('Distance')
xlabel('Covariance')


