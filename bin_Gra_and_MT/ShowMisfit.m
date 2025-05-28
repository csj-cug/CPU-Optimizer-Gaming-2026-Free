figure
    set(gcf,'color','w')  
plot(1:n,misfit1(1:n),'linewidth',2);
hold on
%plot(1:n,misfit2(1:n),'linewidth',2);
set(gca,'fontweight','bold');
xlabel('Iteration number');
xticks(1:1:5)
ylabel('Normalized misfit');
set(gca,'fontweight','bold');
%legend('Gravity data misfit','Magnetic data misfit');