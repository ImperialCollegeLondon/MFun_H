function Z = imagesc_CAT(Z,varargin)
% imagesc Z for a defined Catchment
%
% Input (Mandatory)
% 
% Inputs (Optional Parameter-value pairs)
% 
% Examples
% 
% Author Yuting CHEN
%        Imperial College London
%        yuting.chen17@imperial.ac.uk

p = inputParser;
p.addRequired('Z');
p.addParameter('Cat',[],@ischar);
p.addParameter('RG',[],@ischar);
p.addParameter('caxis',[],@isnumeric);
p.addParameter('axesH',[],@(x)isa(x,'matlab.graphics.axis.Axes'));
p.addParameter('title',[],@ischar);

p.parse(Z,varargin{:});

if ~isempty(p.Results.axesH)
    axes(p.Results.axesH);
else
    figure;
end


if ~isempty(p.Results.Cat)
    load(p.Results.Cat,'catBound','xyMat');
    in = inpolygon(xyMat.xx,xyMat.yy,catBound.x,catBound.y);
    Z(~in) = NaN;
    imAlpha = ones(size(Z));
    imAlpha(isnan(Z)) = 0;
    imagesc(xyMat.x_yr,xyMat.y_yr,Z,'AlphaData',imAlpha);
else
    imAlpha = ones(size(Z));
    imAlpha(isnan(Z)) = 0;
    imagesc(Z,'AlphaData',imAlpha);
end
set(gca,'YDir','normal')
hold on;

if ~isempty(p.Results.Cat)
    plot(catBound.x,catBound.y,'k-');
    hold on;
    plot(catBound.SN.x,catBound.SN.y,'.','markersize',2,'color',[0.2 0.2 0.2]);
    hold on;
end

if ~isempty(p.Results.RG)
    load(p.Results.RG,'rgSite');
    rglh = plot(rgSite.x,rgSite.y,'bo','markersize',3);
    wideM = 5000;
    xlim([min(xyMat.x_yr)-wideM,max(xyMat.x_yr)+wideM]);
    ylim([min(xyMat.y_yr)-wideM,max(xyMat.y_yr)+wideM]);
    legend(rglh,'Rain Gauge','Location','best');
end

axis off
axis tight
axis equal



if ~isempty(p.Results.caxis)
    caxis(p.Results.caxis)
end

if ~isempty(p.Results.title)
    title(p.Results.title)
end


end