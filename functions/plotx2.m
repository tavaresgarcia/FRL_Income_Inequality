function plotx2(t,y)
set(gcf,'DefaultAxesColorOrder',[0 0 0;0 0 0;0 0 0;0 0 0]);
set(gcf,'color','w');
cu=y(:,2);
cl=y(:,3);

h=t;
h=h';
hh=fill([h(1); h(1:end); flipud([h(1:end); h(end)])],[cu(1); cl(1:end); flipud([cu(1:end); cl(size(cl,1))])],'b');
set(hh,'edgecolor',[0.2 0.2 0.2]);
set(hh,'facecolor',[.9 .9 .9]);


hold on
plot(h,y(:,1),'LineWidth',1.3);
axis tight
 hold on;
zz=zeros(size(y,1),1);
 plot(h,zz,'r-');
