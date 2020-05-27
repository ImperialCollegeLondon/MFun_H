classdef RainfallDataClass
    % define RainfallDataClass saving all necessary info for efficient processing
    % @ yuting chen
    % Imperial College London
    %
    %
    % Example:
    % obj = RainfallDataClass(Val,Missingval,ScaleF,Unit,Time,XX,YY,Attributes)
    % 
    % oriDataVal = originalData(obj,'double')
    % 
    % meanVal = nanmean(obj,'spatial')
    % 
    % meanVal = nanmean(obj,'timeseries')
    %
    % war = STATS(obj,'war')
    %
    % subObj = extractOnePeriod(obj,isInPeriod)
    
    properties
        Val
        Missingval
        ScaleF
        Unit % 'mm/h' or ..
        Time % <datatime>
        XX % [m]
        YY % [m]
        Attributes
        dt % [minutes]
        dx % [km]
    end
    methods
        % constructor
        function obj = RainfallDataClass(Val,Missingval,ScaleF,Unit,Time,XX,YY,Attributes)
            arguments
                Val = ones(5,5,10);
                Missingval (1,1) double = NaN;
                ScaleF (1,1) double = 1;
                Unit (1,:) char = 'mm/h';
                Time (:,1) = ones(10,1);
                XX (:,:) double = rand(5,5);
                YY (:,:) double = rand(5,5);
                Attributes (1,:) char = 'none'
            end
            obj.Val = Val;
            obj.Missingval = Missingval;
            obj.ScaleF = ScaleF;
            obj.Unit = Unit;
            obj.Time = Time;
            obj.XX = XX;
            obj.YY = YY;
            obj.Attributes = Attributes;
            obj.dt = minutes(Time(2)-Time(1));
            obj.dx = (YY(2)-YY(1))/1000;
        end
        
        function oriDataVal = originalData(obj,format)
            arguments
                obj
                format char = 'double'
            end
            % return original Data without any formatting/rounding/scaling
            oriDataVal = double(obj.Val);
            oriDataVal(oriDataVal == obj.Missingval) = NaN;
            oriDataVal = oriDataVal/obj.ScaleF; 
        end
        
        function raintime = getTime(obj)
            raintime = obj.Time;
        end
        
        function obj = squeezeData(obj,format)
            arguments
                obj
                format char = 'int16'
            end
            % return original Data without any formatting/rounding/scaling
            obj.Val = obj.Val*obj.ScaleF;
            obj.Val(isnan(obj.Val)) = obj.Missingval;
            obj.Val = int16(obj.Val);
        end
        
        function meanVal = nanmean(obj,dim)
            % return nanmean of obj.Val on 'dim' dimension
            oriDataVal = originalData(obj,'double');
            switch(upper(dim))
                case 'SPATIAL'
                    meanVal = squeeze(nanmean(oriDataVal,3));
                case 'TIMESERIES'
                    meanVal = squeeze(nanmean(reshape(...
                        oriDataVal,[],size(oriDataVal,3)),1));
            end
        end
        
        function stats = STATS(obj,property)
            % calculate specific statsitical property for obj data
            % available property include:
            %            war
            %            imf
            %            dryDur
            %            wetDur
            %            ...
            switch(upper(property))
                case 'WAR'
                    stats = nanmean(reshape(obj.Val,numel(obj.XX),[]) ...
                        > (0.1*obj.ScaleF),1);
                case 'IMF'
                    stats = nanmean(reshape(originalData(obj,'double'),...
                        numel(obj.XX),[]),1);
                case 'CVPOS'
                    % cv in positive part
                    oridata = originalData(obj,'double');
                    oridata = oridata(oridata>0.1);
                    stats = nanstd(oridata(:))./nanmean(oridata(:));
            end
        end
        
        function subObj = extractOnePeriod(obj,isInPeriod)
            % extract some sepecific time period from obj data.
            subObj = obj;
            subObj.Val = squeeze(subObj.Val(:,:,isInPeriod));
            subObj.Time = squeeze(subObj.Time(isInPeriod));
        end
        
        function subObj = extractOneArea(obj,isInArea,XX)
            % extract some sepecific time period from obj data.
            subObj = obj;
            subObj.XX = reshape(subObj.XX(isInArea),size(XX));
            subObj.YY = reshape(subObj.YY(isInArea),size(XX));
            subObj.Val = reshape(subObj.Val,[],size(subObj.Val,3));
            subObj.Val = subObj.Val(isInArea(:),:);
            subObj.Val = reshape(subObj.Val,[size(subObj.XX),size(obj.Val,3)]);
        end
        
    end
end