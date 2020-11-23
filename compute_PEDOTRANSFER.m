function [Osat,L,Pe,Ks,O33] = compute_PEDOTRANSFER(SAND,CLAY,ORG);
% COMPUTE_PEDOTRANFER calculate several soil properties for input of matrix.
% Input: SAND: Matrix<double> [%]
%        CLAY: Matrix<double> [%]
%        ORG:  Matrix<double> [%]
% Output:Osat:
%        L:
%        Pe:
%        Ks:
%        O33:
%
% No-data value: input is set as NaN
%                output is set as NaN
%
% by TitianChen
% Imperial College London

fprintf('---------------Input Statistics--------------\n');
fprintf('Val: SAND       CLAY        ORG\n');
fprintf('Max: %.2f\t\t%.2f\t\t%.2f\n',...
    max(SAND(:)),max(CLAY(:)),max(ORG(:)));
fprintf('Min: %.2f\t\t%.2f\t\t%.2f\n',...
    min(SAND(:)),min(CLAY(:)),min(ORG(:)));
fprintf('Mean %.2f\t\t%.2f\t\t%.2f\n',...
   nanmean(SAND(:)),nanmean(CLAY(:)),nanmean(ORG(:)));
fprintf('...PEDOTRANSFER ing....\n');

[Osat,L,Pe,Ks,O33] = deal(NaN(size(CLAY)));

return

for i = 1:size(CLAY,1)
    for j = 1:size(CLAY,2)
        [Osat(i,j),L(i,j),Pe(i,j),Ks(i,j),O33(i,j),~,~,~,~,~] = ...
            Soil_parameters(SAND(i,j),CLAY(i,j),ORG(i,j));
    end
end
fprintf('...PEDOTRANSFER finished.\n');
fprintf('--------------Summarry Statistics-------------\n');
fprintf('Val: Osat       L           Pe          Ks          O33\n');
fprintf('Max: %.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n',...
    max(Osat(:)),max(L(:)),max(Pe(:)),max(Ks(:)),max(O33(:)));
fprintf('Min: %.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n',...
    min(Osat(:)),min(L(:)),min(Pe(:)),min(Ks(:)),min(O33(:)));
fprintf('Mean %.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n',...
   nanmean(Osat(:)),nanmean(L(:)),nanmean(Pe(:)),nanmean(Ks(:)),nanmean(O33(:)));


end

