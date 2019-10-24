function RES = rigrid_MATs_inStruct(STRUCT,newSize,varargin)

% RIGRID_MATS_INSTRUCT finish the regrid for a struct storing several Matrixes of
% soil properties.
% All field will be maintained and change into the input cellsize
% 
% Input: STRUCT
%        newSize
%        varargin::fillNan <true/false>
% Output:Soil_input_2
% 
% Struct:STRUCT: can have: matrix field, single value, vector.
%             Notice: must have:
%                      STRUCT.xRange: [min,max]
%                      STRUCT.yRange: [min,max]
% 
% by Yuting CHEN

RES = struct;
NAMES = fieldnames(STRUCT);

if ~isempty(NAMES) && iscell(NAMES)
    
    try
        xRange = STRUCT.xRange;
        yRange = STRUCT.yRange;
    catch
        returnNotice;
        return;
    end
    
    for fi = 1:length(NAMES)
        
        MAT = getfield(STRUCT,NAMES{fi}); %#ok<*GFLD>
        CS = 0;
        if size(MAT,2)>2
            MAT_new = imresize(MAT,newSize);
            if ~isempty(varargin) && varargin{1} == true
                MAT_new(isnan(MAT_new)) = nanmean(MAT_new(:));
            end
            % x_aux = linspace(xRange(1),xRange(2),size(MAT,2));
            % y_aux = linspace(yRange(1),yRange(2),size(MAT,1));
            % [xx,yy] = meshgrid(x_aux,y_aux);
            %
            % F = scatteredInterpolant(xx(:),yy(:),MAT(:));
            % Y = linspace(yRange(1),yRange(2),newSize(1));
            % X = linspace(xRange(1),xRange(2),newSize(2));
            % Y = Y(end:-1:1);
            %
            % [xx,yy] = meshgrid(X,Y);
            % MAT_new = F(xx,yy);
            %
            % CS = unique(abs(diff(X)));
            
            RES = setfield(RES,NAMES{fi},MAT_new); %#ok<SFLD>
            
        else
            
            % not matrix --> copy & paste;
            RES = setfield(RES,NAMES{fi},MAT); %#ok<SFLD>
            
        end
    end
    
    fprintf('* RIGRID* Finished.\n Current Cell Size: %.2f.\n',CS);
    
    RES.xRange = [min(xRange),max(xRange)];
    RES.yRange = [min(yRange),max(yRange)];
    
else
    returnNotice;
    return;
end

end

function returnNotice

fprintf('check input STRUCT.\n');

end



