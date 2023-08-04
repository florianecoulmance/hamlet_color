function scatterplot_2D(Y,thesespecies,colors)
% 2D pairwise plots
figure
subplot(131)
gscatter(Y(:,1),Y(:,2),thesespecies',colors,'.',30)
xlabel('Dim1')
ylabel('Dim2')
axis square
set(gca,'fontname','Arial')
set(gca,'fontsize',20)
legend off

subplot(132)
gscatter(Y(:,2),Y(:,3),thesespecies',colors,'.',30)
xlabel('Dim2')
ylabel('Dim3')
axis square
set(gca,'fontname','Arial')
set(gca,'fontsize',20)
legend off

subplot(133)
gscatter(Y(:,1),Y(:,3),thesespecies',colors,'.',30)
xlabel('Dim1')
ylabel('Dim3')
axis square
set(gca,'fontname','Arial')
set(gca,'fontsize',20)
m = get(gca,'legend');
set(m,'location','southeast')
set(m,'fontsize',10)


set(gcf,'units','normalized')
set(gcf,'outerposition',[0 0 1 1])