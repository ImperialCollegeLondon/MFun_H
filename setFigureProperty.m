function handle = setFigureProperty(varargin)
% SETFIGUREPROPERTY(VARARGIN) customize the figure and set the size,
% font...
%
% Input:
%    null
%    nargin 1
%       'Full'    :
%       'High'    :
%       'Wide'    :
%       'Paper'   :
%       'Subplot2':
%       'Subplot3':
%    nargin pairs
%       property
%
%
%
% Examples:
%
%
% by Yuting Chen
%  Imperial College London

set(0,'defaultAxesFontSize',14,'defaultAxesFontName','Cambria','defaultAxesTitleFontWeight','Normal');
if nargin == 0
    set(0,'defaultAxesFontSize',14,'defaultAxesFontName','Cambria','defaultAxesTitleFontWeight','Normal');
    XYWH = [150,150,250,180];
    set(gcf,'units','points','position',XYWH);
else
    if nargin == 1
        if isnumeric(varargin{1})
            set(0,'defaultAxesFontSize',14,'defaultAxesFontName',...
                'Tw cen MT','defaultAxesTitleFontWeight','Normal');
            XYWH = [50,-50,varargin{1}(1),varargin{1}(2)];
            set(gcf,'units','points','position',XYWH);
        else
            switch(varargin{1})
                case 'Big'
                    set(0,'defaultAxesFontSize',14,'defaultAxesFontName','Cambria','defaultAxesTitleFontWeight','Normal');
                    
                    XYWH = [50,-50,400,350];
                    set(gcf,'units','points','position',XYWH);
                case 'Full'
                    set(gcf, 'Position', get(0, 'Screensize'));
                case 'High'
                    XYWH = [50,-50,300,600];
                    set(gcf,'units','points','position',XYWH);
                case 'Wide'
                    XYWH = [50,-50,900,380];
                    set(gcf,'units','points','position',XYWH);
                case 'Paper'
                    set(0,'defaultAxesFontSize',12,'defaultAxesFontName','Cambria','defaultAxesTitleFontWeight','Normal');
                    XYWH = [50,-50,600,450];
                    set(gcf,'units','points','position',XYWH);
                case 'Paper_2'
                    set(0,'defaultAxesFontSize',14,'defaultAxesFontName','Cambria','defaultAxesTitleFontWeight','Normal');
                    XYWH = [50,-50,600,225];
                    set(gcf,'units','points','position',XYWH);
                case 'Paper_4'
                    set(0,'defaultAxesFontSize',14,'defaultAxesFontName','Cambria','defaultAxesTitleFontWeight','Normal');
                    XYWH = [50,-50,600,450];
                    set(gcf,'units','points','position',XYWH);
                    
                case 'Subplot4'
                    set(0,'defaultAxesFontSize',14,'defaultAxesFontName','Cambria','defaultAxesTitleFontWeight','Normal');
                    XYWH = [50,-50,750,225];
                    set(gcf,'units','points','position',XYWH);
                    
                case 'Subplot2'
                    set(0,'defaultAxesFontSize',14,'defaultAxesFontName','Cambria','defaultAxesTitleFontWeight','Normal');
                    XYWH = [50,-50,550,220];
                    set(gcf,'units','points','position',XYWH);
                case 'Subplot3'
                    set(0,'defaultAxesFontSize',14,'defaultAxesFontName','Cambria','defaultAxesTitleFontWeight','Normal');
                    XYWH = [50,-50,750,220];
                    set(gcf,'units','points','position',XYWH);
                otherwise
            end
        end
    else
        
    end
end
handle = gcf;
end