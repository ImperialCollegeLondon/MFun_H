%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Input:
%sand: sand fraction [%], (nx1)
%clay: Clay fraction[%], (nx1)
%T: type of classification: 'usda'
%Varargin{1}: PLOT, BOOLEAN: 1: plots the usda texture triangle plus data
%points
%
%Output:
%CellArray of soil classes
%silt: silt fraction [%], (nx1)
%
%Example:
%soil_classification(sand,clay, 'usda')
%soil_classification([0.5;0.3;0.1;0],[0.1;0.2;0.4;0.8], 'usda',1)
%--------------------------------------------------------------------------
%disp('-------------------------------------------------------------------')
%disp('                       U S D A soil Classification                 ')
%disp('-------------------------------------------------------------------')
function [SC, silt]=soil_classification(sand, clay, T,varargin)
 if strcmp(T,'usda')==0
    error('only usda defined so far...')
 end
 if isempty(sand)|isnan(sand)|min(sand)<0
     error('wrong input: sand')
 end
 if isempty(clay)|isnan(clay)|min(clay)<0
     error('wrong input: sand')
 end
 if max(sand)>1
     sand=sand/100;
     warning('sand is >1; ...divided by 100 to get fraction instead of %')
 end
 if max(clay)>1
     clay=clay/100;
     warning('clay is >1; ...divided by 100 to get fraction instead of %')
 end
 if max(sand+clay)>1
    error('data inconsisent: (Clay + Sand)>1')
 end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
silt=1-sand-clay;
i=1;
 for i=1:length(clay)
     if (silt(i)+1.5*clay(i))<=0.15
         SC{i,:}='SAND';
     elseif ((silt(i)+1.5*clay(i))>=.15)&((silt(i)+2*clay(i))<=0.3)
         SC{i,:}='LOAMY SAND';
     elseif (clay(i)>=0.07) & (clay(i)<=0.2) & (sand(i)>=0.52) & ((silt(i)+2*clay(i))>=0.3) 
         SC{i,:}='SANDY LOAM';
     elseif (clay(i)<=0.07) & (silt(i) <= 0.5) & ((silt(i)+2*clay(i))>=0.3)
         SC{i,:}='SANDY LOAM';
     elseif (clay(i)>=0.07) & (clay(i)<=0.27) & (silt(i)>=0.28) & (silt(i)<=0.5) & (sand(i)<=0.52)
         SC{i,:}='LOAM';
     elseif ((silt(i)>=0.5) & clay(i)>0.12 & clay(i)<0.27) | (silt(i)>=0.5 & silt(i)<=0.8 & clay(i)<=0.12)
         SC{i,:}='SILT LOAM';
     elseif silt(i)>=0.8 & clay(i)<=0.12
         SC{i,:}='SILT';
     elseif clay(i)>=0.2 & clay(i)<=0.35 & silt(i)<=0.28 & sand(i)>=0.45
         SC{i,:}='SANDY CLAY LOAM';
     elseif clay(i)>=0.27 & clay(i)<=0.4 & sand(i)>=0.2 & sand(i)<=0.45
         SC{i,:}='CLAY LOAM';
     elseif clay(i)>=0.27 & clay(i)<=0.4 & sand(i)<=0.2
         SC{i,:}='SILTY CLAY LOAM';
     elseif clay(i)>=0.35 & sand(i)>=0.45
         SC{i,:}='SANDY CLAY';
     elseif clay(i)>=0.4 & silt(i)>=0.4
         SC{i,:}='SILTY CLAY';
     elseif clay(i)>= 0.4 & sand(i)<=0.45 & silt(i)<=0.4
         SC{i,:}='CLAY';
     else
         error('CHECK!Function Soil_Classification.m')
     end
 end
 %------------------------------------------------------------------------
 if varargin{1}==1
    %close
     %set(gcf,'Color','w','position',[2800 300 600 600])
     plot([0 100],[0 0],'k');hold on
     set(gca,'XDir','reverse')
     plot([0 50],[0 100],'k')
     plot([50 100],[100 0],'k')
     xlabel('Sand [%]')
     text(95,60,'Clay [%]','FontSize',12)
     text(25,60,'Silt [%]','FontSize',12)
     text(55,-5,'Sand [%]','FontSize',12)

     p1=fill([100 85 95 100],[0 0 10 0],[0.3 0.3 0.3]);
     text(97,3, 'Sand','color','w')
     
     p2= fill([85 70 92.4 95 85],[0 0 15 10 0],'y');
     text(87,3,{'Loamy';'   sand'})
     
     p3=fill([70 50 46.7 55 62 90 92.4 70],[0 0 7 7 20 20 15 0],[0.6 0.3 0]);
     text(75,10,'Sandy Loam')    
     
     p4=fill([55 46.7 38 58.1 62 55],[7 7 27 27 20 7],[0 0.7 0]);
     text(55,15,'Loam')
         
     p5=fill([50 20 13.5 6 13.2 38 50],[0 0 12 12 27 27 0],[0.5 0.5 0]);
     text(35, 10, 'Silt Loam')
     p6=fill([20 0 6 13.5 20],[0 0 12 12 0],[0.3 0.6 0]);
     text(12,5,'Silt')
     
     p7=fill([90 62 58.1 62.5 82.5 90],[20 20 27 35 35 20],[0.9 0.3 0.1]);
     text(85,27,'Sandy Clay Loam')

     p8=fill([58.1 33 40 65 58],[27 27 40 40 27],[0 0.3 0.8]);
     text(55,34,'Clay Loam','Color','w')

     p9=fill([33 13.5 20 40 33],[27 27 40 40 27],[0.8 0.8 0.8]);
     text(33,33,{'Silty';'Clay Loam'})
     
     p10=fill([82.5 62.5 72.5 82.5],[35 35 55 35],'k');
     text(80,38,'Sandy Clay','Color','w')
     
     p11=fill([65 40 30 50 72.5 65],[40 40 60 100 55 40],[1 0.8 0.2]);
     text(55,60,'Clay')
     
     p12=fill([40 20 30 40],[40 40 60 40],[0.2 0.2 0.2]);
     text(35,45,'Silty Clay','Color','w')
     axis square
     axis off
     
     o=100-(1-clay)*100; %offset
     plot(sand*100+o/2,clay*100,'ok','MarkerFaceColor','w','Markersize',3)
 end
end %of function