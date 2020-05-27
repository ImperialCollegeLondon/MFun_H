function format_xylabel(ha,rowno,colno)
set(ha(((1:rowno)-1)'*colno+(2:colno)),'YTickLabel',[],'YLabel',[])
set(ha(1:(colno*(rowno-1))),'XTickLabel',[],'XLabel',[])
end