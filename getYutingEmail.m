function str = getYutingEmail(varargin)

if isempty(varargin)
    str = 'yutingchen0604@hotmail.com';
    
else
    switch(upper(varargin{1}))
        case('PASSWORD')
            str = 'AaBb14207';
        otherwise
            
    end
end

end