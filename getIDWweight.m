function weight = getIDWweight(DIST,P)
%
% Inverse Distance Weighting method

dis = DIST(:);
weight = ((1./(dis.^P))./nansum(1./(dis.^P)));
weight = reshape(weight,size(DIST));

end