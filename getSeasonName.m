function [str] = getSeasonName(season)
switch(season)
    case 1
        str = 'MAM';
    case 2
        str = 'JJA';
    case 3
        str = 'SON';
    case 4
        str = 'DJF';
    otherwise
end
end

