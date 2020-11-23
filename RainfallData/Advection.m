classdef Advection
    % This class is to define advection information of rainfall images
    % @ Yuting Chen
    
    properties
        Vx % [km/h]
        Vy % [km/h]
        Vstd % [km/h] (over space).
        Vdir % [-pi~pi]
        unit = 'Km/h'
    end
    methods
        function [obj] = smooth(obj,smoothN)
            obj.Vx = smooth(obj.Vx,smoothN);
            obj.Vy = smooth(obj.Vy,smoothN);
        end
        function [obj] = updateAdvection(obj,Vx,Vy,Vstd,Vdir)
            obj.Vx = Vx;
            obj.Vy = Vy;
            obj.Vstd = Vstd;
            obj.Vdir = Vdir;
        end
        function [obj] = Advection(obsDATA,filter)
            arguments
                obsDATA (1,1) RainfallDataClass = RainfallDataClass()
                filter (1,1) double = 0.1;
            end
            
            RE = originalData(obsDATA,'double');
            
            VE = compute4OneStorm(RE);
            
            obj.Vx = VE.rspeedx;
            obj.Vy = VE.rspeedy;
            obj.Vstd = VE.rspeedstd;
            obj.Vdir = VE.rspeeddir;
            
            function stats = compute4OneStorm(Rain)
                % unit:  Km/h
                dt = obsDATA.dt/60;
                dx = obsDATA.dx;
                
                [rspeedx,rspeedy,rspeedstd,rspeeddir] = deal(NaN);
                
                zeroVal = 0;
                if size(Rain,3)==length(obsDATA.Time)
                    Rtemp = permute(Rain,[3,1,2]);% change 'time steps' into the first dimension;
                end
                Rtemp = func_filter(Rtemp,filter,zeroVal);%%%%%%%%%
                
                opticFlow = opticalFlowLK('NoiseThreshold',0.0005);
                
                R0 = squeeze(Rtemp(1,:,:));
                convDn = 120*dt;% corresponds to maximum movement within this timestep
                R0 = conv2(R0,ones(convDn)/convDn/convDn,'same');
                R0 = func_R2dBZ(R0,'UKMO',zeroVal);
                for sni = 1:size(Rtemp,1)
                    
                    R1 = squeeze(Rtemp(sni,:,:));
                    R1 = conv2(R1,ones(convDn)/convDn/convDn,'same');
                    R1 = func_R2dBZ(R1,'UKMO',zeroVal);
                    
                    flow = estimateFlow(opticFlow,R1);
                    
                    V = dx/dt * nanmean(flow.Magnitude(R1~=zeroVal | R0~=zeroVal));
                    
                    Vx = dx/dt * nanmean(flow.Vx(R1~=zeroVal | R0~=zeroVal));
                    Vy = dx/dt * nanmean(flow.Vy(R1~=zeroVal | R0~=zeroVal));
                    Vstd = dx/dt * nanstd(flow.Magnitude(R1~=zeroVal | R0~=zeroVal));
                    Vdir = nanstd(flow.Orientation(R1~=zeroVal | R0~=zeroVal));
                    
                    if sni > 1
                        rspeedx(sni) = Vx;
                        rspeedy(sni) = Vy;
                        rspeedstd(sni) = Vstd;
                        rspeeddir(sni) = Vdir;
                    end
                    
                    %             if (~any(R0(:)~=zeroVal) & any(R1(:)~=zeroVal))|...
                    %                     (~any(R1(:)~=zeroVal) & any(R0(:)~=zeroVal))%#ok<OR2,AND2>
                    %                 % In the case that: having peak(>5mm/h) pixels in image 0/1 but no peak
                    %                 % in image 1/0, the speed is unknown (because the things
                    %                 % happen within one hour which is not able to be captured
                    %                 % in current output.
                    %                 rspeedx(sni) = NaN;
                    %                 rspeedy(sni) = NaN;
                    %             end
                    % hold off;
                    R0 = R1;
                end
                
                rspeedx = rspeedx(:);
                rspeedy = rspeedy(:);
                rspeedstd = rspeedstd(:);
                rspeeddir = rspeeddir(:);
                
                stats = table(rspeedx,rspeedy,rspeedstd,rspeeddir);
            
                % plot(rspeed);
                % drawnow
                % pause(0.2);
                
                function Rtemp = func_filter(Rtemp,thres,zeroVal);
                    Rtemp(Rtemp<thres) = zeroVal;
                end
                function Rtemp = func_R2dBZ(Rtemp,source,zeroVal)
                    if strcmp(source,'UKMO')
                        Rtemp(Rtemp~=zeroVal) = (log10(200)+1.6*log10(Rtemp(Rtemp~=zeroVal)))*10;
                    end
                end
            end
        end
    end
end