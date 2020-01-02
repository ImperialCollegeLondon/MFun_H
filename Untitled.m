% save PET and GEAR rain for one location
X = 282900;
Y = 283800;
YEAR = 2000:2010;

%%
[PET,status] = importCHESS_PET(XX,YY,YEAR);
PET = squeeze(PET);

%%
t = datetime(2000,1,1):1:datetime(2010,12,31);
t = t';
PR = [];
%%
for year = 2010%YEAR
date0 = datetime(year,1,1):1:datetime(year,12,31);
[~,~,RAIN] = import_GEAR_DAILY(X,Y,datenum(date0),0);
PR = [PR;reshape(squeeze(RAIN),[],1)];
end

%%
save('data.mat','PET','PR','X','Y','t');
