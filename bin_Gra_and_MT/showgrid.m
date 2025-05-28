figure;
hp=pcolor(squeeze(X(:,2,:)),squeeze(Z(:,2,:)),0*squeeze(X(:,2,:)));
xlabel('x');
ylabel('z');
set(gca,'ydir','reverse');
hold on
% fill([xa1 xa1 xa2 xa2 xa1],[za1 za2 za2 za1 za1],'b')
%fill([min(xa{1}) min(xa{1}) max(xa{1}) max(xa{1}) min(xa{1})],[min(za{1}) max(za{1}) max(za{1}) min(za{1}) min(za{1})],'b')
%fill([min(xa{2}) min(xa{2}) max(xa{2}) max(xa{2}) min(xa{2})],[min(za{2}) max(za{2}) max(za{2}) min(za{2}) min(za{2})],'b')
%fill([min(xa{3}) min(xa{3}) max(xa{3}) max(xa{3}) min(xa{3})],[min(za{3}) max(za{3}) max(za{3}) min(za{3}) min(za{3})],'b')
%fill([min(xa{4}) min(xa{4}) max(xa{4}) max(xa{4}) min(xa{4})],[min(za{4}) max(za{4}) max(za{4}) min(za{4}) min(za{4})],'b')
set(gca,'fontweight','bold');
xlabel('x(m)');
ylabel('Z(m)');
hold on
% plot(xs,zs,'r*','markersize',10);
set(hp,'EdgeAlpha',.5);

figure;
hp=pcolor(squeeze(Y(2,:,:)),squeeze(Z(2,:,:)),0*squeeze(Y(2,:,:)));
xlabel('Y');
ylabel('z');
set(gca,'ydir','reverse');
hold on;
% fill([ya1 ya1 ya2 ya2 ya1],[za1 za2 za2 za1 za1],'b')
%fill([min(ya{1}) min(ya{1}) max(ya{1}) max(ya{1}) min(ya{1})],[min(za{1}) max(za{1}) max(za{1}) min(za{1}) min(za{1})],'b')
%fill([min(ya{2}) min(ya{2}) max(ya{2}) max(ya{2}) min(ya{2})],[min(za{2}) max(za{2}) max(za{2}) min(za{2}) min(za{2})],'b')
%fill([min(ya{3}) min(ya{3}) max(ya{3}) max(ya{3}) min(ya{3})],[min(za{3}) max(za{3}) max(za{3}) min(za{3}) min(za{3})],'b')
%fill([min(ya{4}) min(ya{4}) max(ya{4}) max(ya{4}) min(ya{4})],[min(za{4}) max(za{4}) max(za{4}) min(za{4}) min(za{4})],'b')
set(gca,'fontweight','bold');
xlabel('Y(m)');
ylabel('Z(m)');
hold on
% plot(ys,zs,'r*','markersize',10);
set(hp,'EdgeAlpha',.5);


figure;
hp=pcolor(squeeze(X(:,:,2)),squeeze(Y(:,:,2)),0*squeeze(X(:,:,2)));
xlabel('x');
ylabel('y');
hold on
% fill([xa1 xa1 xa2 xa2 xa1],[ya1 ya2 ya2 ya1 ya1],'b')
%fill([min(xa{1}) min(xa{1}) max(xa{1}) max(xa{1}) min(xa{1})],[min(ya{1}) max(ya{1}) max(ya{1}) min(ya{1}) min(ya{1})],'b')
%fill([min(xa{2}) min(xa{2}) max(xa{2}) max(xa{2}) min(xa{2})],[min(ya{2}) max(ya{2}) max(ya{2}) min(ya{2}) min(ya{2})],'b')
%fill([min(xa{3}) min(xa{3}) max(xa{3}) max(xa{3}) min(xa{3})],[min(ya{3}) max(ya{3}) max(ya{3}) min(ya{3}) min(ya{3})],'b')
%fill([min(xa{4}) min(xa{4}) max(xa{4}) max(xa{4}) min(xa{4})],[min(ya{4}) max(ya{4}) max(ya{4}) min(ya{4}) min(ya{4})],'b')
set(gca,'fontweight','bold');
xlabel('X(m)');
ylabel('Y(m)');
hold on
% plot(xs,ys,'r*','markersize',10);
set(hp,'EdgeAlpha',.5);


figure;
m1  = reshape(sigx,60,60,46);
SZ1 = reshape(m1(:,:,15),[60,60]).';
s1  = imagesc(SZ1);colormap(flipud(jet));
 h  = colorbar('south');
