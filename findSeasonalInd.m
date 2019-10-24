function [season_ind] = findSeasonalInd(month)
season_ind = cell(1,4);
season_ind{1} = find(month==3 | month==4 | month==5);% MAM
season_ind{2} = find(month==6 | month==7 | month==8);% JJA
season_ind{3} = find(month==9 | month==10 | month==11);% SON
season_ind{4} = find(month==12 | month==1 | month==2);% DJF
end