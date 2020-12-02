%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SPEARMAN'S RHO TEST   %%%%% %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Td,p_value]=SpearmanRho(V,alpha)
%%%%%%%%%%%%%%%%%
%%% Performs Spearman's rho test of the null hypothesis of trend
%%% absence in the vector V,  against the alternative of trend. 
%%% The result of the test is returned in Td = 1 indicates positive trend
%%% Td = -1 indicates negative trends, i.e., 
%%% a rejection of the null hypothesis at the alpha significance level. Td = 0 indicates
%%% a failure to reject the null hypothesis at the alpha significance level.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUTS
%V = time series [vector]
%alpha =  significance level of the test [scalar]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% From Matlab Help %%%%%%%%%%%%%%%
%The significance level of a test is a threshold of probability a agreed
%to before the test is conducted. A typical value of alpha is 0.05. If the p-value of a test is less than alpha,
%the test rejects the null hypothesis. If the p-value is greater than alpha, there is insufficient evidence 
%to reject the null hypothesis. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUTS
%Td = test result: [1] Positive Trend; [-1] Negative Trend;  [0] Insufficient evidence to reject the null hypothesis
%p_value = p-value of the test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% From Matlab Help %%%%%%%%%%%%%%%
%The p-value of a test is the probability, under the null hypothesis, of obtaining a value
%of the test statistic as extreme or more extreme than the value computed from
%the sample.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% References 
%Daniels, H. E. (1950), Rank correlation and population models, J. R. Stat.
% Soc. Ser. B, 12, 171? 181.
%Khaliq, M. N., T. B. M. J. Ouarda, P. Gachon, L. Sushama, and A. St-Hilaire
%(2009), Identification of hydrological trends in the presence of serial and 
%cross correlations: A review of selected methods and their application to 
%annual flow regimes of Canadian rivers, J. Hydrol., 368, 117 ? 130, 
%doi:10.1016/j.jhydrol.2009.01.035
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simone Fatichi -- simonef@dicea.unifi.it
%   Copyright 2009
%   $Date: 2009/10/03 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V=reshape(V,length(V),1); 
alpha = alpha/2; %
n=length(V); 
rg_val=1:n;
[V,I]=sort(V);
d=abs(rg_val'-I); 
ro=1-6*sum(d.^2)/(n^3-n); 
if ro == 1 
    tcamp = 0; 
else
tcamp=ro*sqrt((n-2)/(1-ro^2));  %% n-2 dof 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tcrit=tinv(1-alpha,n-2);  
H1=abs(tcamp) > tcrit; %%
Td=sign(tcamp)*H1; 
p_value=2*(1-tcdf(abs(tcamp),n-2)); 
return 