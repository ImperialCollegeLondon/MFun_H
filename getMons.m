function [monv] = getMons(season)
%
% 2:JJA
%

switch(season)
    case 4
        monv = [12 1 2];
    case 1
        monv = [3 4 5];
    case 2
        monv = [6 7 8];% JJA
    case 3
        monv = [9 10 11];
end

end