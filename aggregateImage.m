function [OUT] = aggregateImage(IN,num,option)

arguments
    IN
    num
    option (1,:) char = 'mean'

end
if num == 1
    OUT = IN;
else
    if strcmpi(option,'max')
        blocksize = [num,num];
        OUT = blockproc(IN, blocksize, @(x)max(x.data(:)));
    else
        OUT = imresize(IN,1/num,'method','box');
    end
end

end