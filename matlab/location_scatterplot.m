function location_scatterplot(Y,theselocs,titleWord)


gscatter(Y(:,1),Y(:,2),theselocs',[],'.',30)
m = get(gca,'legend');
set(m,'location','southeast')

axis square
set(gca,'fontsize',40)
set(gca,'fontname','Arial')
xlabel('dim1')
ylabel('dim2')
set(gcf,'units','normalized')
set(gcf,'outerposition',[0 0 1 1])
title(titleWord)
