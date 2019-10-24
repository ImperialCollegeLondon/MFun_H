function res = extractPartMatrix(B,I,J)
% EXTRACTSOIL return the struct storing the data from B with a specific range(I,J)
%
% Input: B: struct storing several fields, each of which is a double matrix
%           with the constant matrix size.
%        I: vector including starting-ending point
%        J: vector inclusing starting-ending point
% Output:res: same structure as B
% 
% 
% struct: (used in B,res)
%        B.*field* <- Matrix<double>
% 
% Example:
%        B = struct;
%        B.a = rand(100,200);
%        B.b = rand(100,200);
%        I = [10,20];
%        J = [20,30];
%        res = extractPartMatrix(B,I,J);
%
% @ Yuting CHEN

res = struct;
NAMES = fieldnames(B);
if ~isempty(NAMES) && iscell(NAMES)
    for fi = 1:length(NAMES)
        F = getfield(B,NAMES{fi});
        mat = F(I(1):I(2),J(1):J(2));
        res = setfield(res,NAMES{fi},mat);
    end
else
    return
end


end




