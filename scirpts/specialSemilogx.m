
function ha =  specialSemilogx(xx,yy)
YLIM = ([min(yy),max(yy)]);
try

ha = tight_subplot(1,2,[0 0],[.1 .01],[.1 .1]);

axes(ha(1))
% left 50
thre = 50;
x1 = xx(xx>thre);
y1 = yy(xx>thre);
hold on;
semilogx(x1,y1,'Linewidth',2);
xlim([50 100]);
set(gca,'XDir','reverse');
ylim(YLIM);

ax.XMinorTick = 'on';

axes(ha(2))
x1 = xx(xx<thre);
y1 = yy(xx<thre);
hold on;
semilogx(x1,y1,'Linewidth',2);
ax = gca;
ax.XAxis.Scale = 'log';
freq = [50 40 30 20 10];
ax.XTick = 50-freq;
ax.XTickLabel = num2cell([freq]);
ax.XMinorTick = 'on';

% ax.XTickLabelMode = 'auto';
set(gca,'XDir','reverse');
ax.YColor = 'none';
xlim([0 50]);
ylim(YLIM);
% xlim([0 50])
catch
    1;
end
end
