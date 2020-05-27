%% This function returns the spatial covariance of the signal val with 
% locations X and Y
% The covariance is estimated at binned intervals with width ds out to a
% maximum distance smax
% The outputs of the function are the spatial covariance values and the
% covariance distance.
function [covariance,distance,r,disI,disJ,Ind_is,Ind_js]=spatialcov(X,Y,D,ds,smax,r,disI,disJ,Ind_is,Ind_js)
%
% Note:
% This script is designed for QUICKLY computing spatial autocovariance from a
% huge amount of images.
% @ Yuting
%

if length(size(D)) == 2
    D = reshape(D,[size(D),1]);
end
if prod(size(D,[1,2]))>12100 || prod(size(D,[1,2,3]))>12100*20000
    error('Normal PCs Memory might be not enough to compute this BIG array in this fast-script.')
end

X = X(:);
Y = Y(:);
val = D(:);
vest=std(val);
disp(['Mean: ',num2str(mean(val)),' Variance: ',num2str((vest))])
distance = NaN(round(smax/ds/2),size(D,3));
covariance = zeros(round(smax/ds)+2,size(D,3));
ncov = zeros(round(smax/ds)+2,1);

if isempty(r)
    r = zeros(length(X),length(X));
    [disI,disJ] = deal(cell(1+smax,1));
    for i = 1:length(X)
        for j = 1:length(X)
            if i <= j % prevent double counting
                r(i,j) = sqrt((X(i)-X(j)).^2+(Y(i)-Y(j)).^2);
            end
        end
    end
    r = int16(round(r));
end
distance = [1:smax]';
[Ind_is,Ind_js] =  deal([]);
for ir = 1:length(distance)
    [is,js] = find(r==distance(ir));
    Ind_is{ir,1} = int16(is);
    Ind_js{ir,1} = int16(js);
end
D = reshape(D,[],size(D,3));
for ir = 1:length(distance)
    tic
    is = Ind_is{ir};
    js = Ind_js{ir};
    thisWholeLen = length(is);
    setStart = 1;
    setLen = 500000;% this number is set based on yt PC's total RAM.
    for setNo = 1:ceil(thisWholeLen/setLen)
        setInd = setStart:1:min(setStart+setLen-1,thisWholeLen);
        if ~isempty(setInd)
            covariance(ir,:) = covariance(ir,:)+sum(D(is,:) & D(js,:),1) ;
            setStart = setInd(end)+1;
        end
    end
    ncov(ir) = ncov(ir)+length(is);
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

for imageNo = 1:size(D,2)
    covariance(:,imageNo) = covariance(:,imageNo)./ncov;  
end
covariance((length(distance)+1):end,:) = [];

% figure
plot(distance,nanmean(covariance,2),'o')
xlabel('Distance')
xlabel('Covariance')




