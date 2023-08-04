function scatterplot_3D(Y,allspecies,uniques,colors)

figure
n = 1;
for i = 1:numel(uniques)
    ind = find(ismember(allspecies,uniques{i}));
    if ~isempty(ind)
    p(n) = plot3(Y(ind,1),Y(ind,2),Y(ind,3),'marker','o','markerfacecolor',colors(i,:),...
        'markeredgecolor','k','color','none','markersize',10);
    hold on
    
    legendvals{n} = [uniques{i},'(',num2str(numel(ind)),')'];
    n = n+1;
    else
    end
end
legend(p,legendvals)

axis square
set(gca,'fontsize',40)
set(gca,'fontname','Arial')
xlabel('MDS Dim 1')
ylabel('MDS Dim 2')
zlabel('MDS Dim 3')
set(gcf,'units','normalized')
set(gcf,'outerposition',[0 0 1 1])


