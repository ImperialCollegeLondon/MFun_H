function str = getMonthName(Mon,outputstr,outputLength)
arguments
    Mon (:,:) %double
    outputstr (1,1) logical = 1
    outputLength (1,1) double = 3;
end

if isempty(Mon)
    str = '';
    return
end

t = datetime(2014,Mon,2);
str = month(t,'shortname');


if outputLength ~= 3
    str = cellfun(@(x)x(1:outputLength),str');
end

if outputstr == true && iscell(str) && length(str) == 1
    str = str{1};
end
end