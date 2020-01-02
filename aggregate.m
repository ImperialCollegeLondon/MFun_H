function out=aggregate(data,scale,method)

if scale==1
    out=data;
else
    A=buffer(data,scale);
    if ischar(method)
        
        if (strcmp(method,'mean')) 
%             out=mean(A);
            out=nanmean(A);
        elseif(strcmp(method,'sum')) 
%             out=sum(A);
            out=nansum(A);
        else
            error('function out=aggregate() BUG 1 !!!!')
        end
        
    else
        
        if  method == 1
%             out=mean(A);
            out=nanmean(A);
        elseif method == 2
%             out=sum(A);
            out=nansum(A);
        else
            error('function out=aggregate() BUG 2 !!!!')
        end
        
    end
    
end
% disp(length(out));
out=out(:);
end
