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
            if strcmp(format,'double')
                oriDataVal = double(obj.Val);
            else
                oriDataVal = single(obj.Val);
            end
            oriDataVal(oriDataVal == obj.Missingval) = NaN;
            oriDataVal = oriDataVal/obj.ScaleF; 
            oriDataVal(oriDataVal < 0) = NaN;
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
            oriDataVal = originalData(obj,'single');
            switch(upper(dim))
                case 'SPATIAL'
                    meanVal = squeeze(nanmean(oriDataVal,3));
                case 'TIMESERIES'
                    meanVal = squeeze(nanmean(reshape(...
                        oriDataVal,[],size(oriDataVal,3)),1));
                case 'MONTH'
                    meanVal = [];
                    for mon = 1:12
                        thisObj = extractOnePeriod(obj,obj.Time.Month == mon);
                        meanVal(mon,1) = STATS(thisObj,'mean');
                    end
                case 'MONTHCV'
                    meanVal = [];
                    for mon = 1:12
                        thisObj = extractOnePeriod(obj,obj.Time.Month == mon);
                        meanVal(mon,1) = STATS(thisObj,'cv');
                    end
                case 'MONTHCVPOS'
                    meanVal = [];
                    for mon = 1:12
                        thisObj = extractOnePeriod(obj,obj.Time.Month == mon);
                        meanVal(mon,1) = STATS(thisObj,'cvpos');
                    end
            end
            meanVal = double(meanVal);
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
                    A = reshape(obj.Val,numel(obj.XX),[]);
                    stats = nansum(A > (0.1*obj.ScaleF), 1)./nansum(A >= 0, 1);
                case 'IMF'
                    stats = double(nanmean(reshape(originalData(obj,'single'),...
                        numel(obj.XX),[]),1));
                case 'MEAN'
                    stats = nanmean(reshape(originalData(obj,'single'),[],1));
                case 'CVPOS'
                    % cv in positive part of distribution of spatial field
                    % for streap computation
                    oridata = originalData(obj,'single');
                    thisStats = [];
                    for ti = 1:size(oridata,3)
                        thisdata = squeeze(oridata(:,:,ti));
                        thisdata = thisdata(thisdata>0.01);
                        if length(thisdata)>10
                            thisStats(ti) = nanstd(thisdata(:))./nanmean(thisdata(:));
                        else
                            thisStats(ti) = NaN;
                        end
                    end
                    stats = double(nanmean(thisStats));
                case 'CV'
                    % cv % for streap computation
                    oridata = originalData(obj,'single');
                    stats = double(nanstd(oridata(:))./nanmean(oridata(:)));
            end
        end
        
        function obj = aggregateDC(obj,scaleDx,scaleDt,Unit)
            if scaleDx ~= 1 || scaleDt ~= 1
                oridata = originalData(obj,'single');
                if logical(~mod(scaleDx,1)) && logical(~mod(scaleDt,1))
                    oridata = imresize3(oridata,'scale',[1/scaleDx,1/scaleDx,1/scaleDt],...
                        'method','box');
                else
                    % #uncomplete#
                    % not 100% precise
                    oridata = imresize3(oridata,'scale',[1/scaleDx,1/scaleDx,1/scaleDt],...
                        'method','box');
                end
                if strcmp(Unit,'mm/day') && strcmp(obj.Unit,'mm/h')
                    oridata = oridata*12;
                end
                obj.Val = int16(oridata*32);
                obj.Val(isnan(oridata)) = obj.Missingval;
                obj.Time = obj.Time(1:scaleDt:end);
                obj.XX = obj.XX(1:scaleDx:end,1:scaleDx:end);
                obj.YY = obj.YY(1:scaleDx:end,1:scaleDx:end);
                obj.dt = obj.dt*scaleDt;
                obj.dx = obj.dt*scaleDx;
            end
        end
        
        function objs = appendTime(obj1,obj2)
            if isempty(obj1)
                objs = obj2;
            else
                objs = obj1;
                objs.Val = cat(3,obj1.Val,obj2.Val);
                objs.Time = cat(1,obj1.Time,obj2.Time);
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
            subObj.Val = reshape(subObj.Val,[size(subObj.XX),size(obj.Val,2)]);
        end
        
    end
end