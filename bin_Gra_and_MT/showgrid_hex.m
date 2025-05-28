function showgrid_hex(X,Z_new,Ztop_interp,xs,ys,zs,xa1,xa2,ya1,ya2,za1,za2)
%% For example: showgrid_hex(X,MeshTop.Z_new,MeshTop.Ztop_interp,-1500,0,900,-500,500,-500,500,2000,2100);
%Z_new and Ztop_interp is stored in MeshTop;
% xa1 xa2 indicates the range of anomalous domain;similar as ya1 ya2 za1
% za2
x = unique(X);
nx = size(X,1);
ind = round(nx/2);

figure;
hp=pcolor(squeeze(X(:,ind,:)),squeeze(Z_new(:,ind,:)),0*squeeze(X(:,ind,:)));
set(gca,'ydir','reverse');
% axis([-2000 2000 00 3000])


xploy = unique(X);
zploy = Ztop_interp(:,ind);
xploy = [xploy;xploy(1)];
zploy = [zploy;zploy(1)];
hold on;fill(xploy,zploy,'b')

hold on;
plot(unique(X),Ztop_interp(:,ind),'r');
hold on;
plot(unique(X),1000*ones(size(x)),'r');

hold on
fill([xa1 xa1 xa2 xa2 xa1],[za1 za2 za2 za1 za1],'b')
set(gca,'fontweight','bold');
xlabel('x(m)');
ylabel('z(m)');
hold on
plot(xs,zs,'r*','markersize',10);
set(hp,'EdgeAlpha',.5);



% C=zeros(87,58,3);
% C(:,:,2)=1;
% figure;surf(squeeze(X(:,1,:)),squeeze(Y(:,1,:)),squeeze(Z(:,1,:)),C);