function ShowInvModelSlice(Sx,Sy,Sz,m1,Nmx,Nmy,Nmz,dmx,dmy,dmz,mx,my,mz)
load D:\Gauss_Newton_Inversion_for_Seismic\gauss_newton_inversion40×40×26\CaseK1\SepInv4\vgrids_true.mat
load vgrids_true vgrids_true
m1         = reshape(m1,Nmx,Nmy,Nmz);
Realm1     = reshape(vgrids_true,Nmx,Nmy,Nmz);
Realm1     = flip(Realm1,3);
m1         = flip(m1,3);
f1         = figure;
dmx        = dmx/1000;
dmy        = dmy/1000;
dmz        = dmz/1000;
mx         = mx/1000;
my         = my/1000;
mz         = mz/1000;
set(f1,'position',[400 311 730 582])


for n = 1:length(Sy)
    subplot(221)
    set(gcf,'color','w')
    SZ1 = reshape(Realm1(:,:,Sz(n)),[Nmx,Nmy]).';
    s1  = imagesc(mx,my,SZ1);colormap(flipud(jet));
    xlabel('X(km)','fontweight','bold','fontsize',15);
    ylabel('Y(km)','fontweight','bold','fontsize',15);
    set(gca,'fontweight','bold','ydir','normal', 'position',[0.1 0.58 0.4 0.4]);
    axis equal;
    axis([min(mx)-dmx/2 max(mx)+dmx/2 min(my)-dmy/2 max(my)+dmy/2])
    title('(a)','position',[-25,75],'fontweight','bold','fontSize',16);
    colorbar('off')
    caxis([2.7,5.8]);
    hold on
    
    subplot(223)
    SY1 = reshape(Realm1(:,Sy(n),:),[Nmx,Nmz]).';
    s2  = imagesc(mx,mz,SY1);colormap(flipud(jet));
    set(gca,'position',[0.1,0.3,0.4,0.2])
    xlabel('X(km)','fontweight','bold','fontsize',15);
    ylabel('Z(km)','fontweight','bold','fontsize',15);
    set(gca,'fontweight','bold', 'ydir','reverse');
    colorbar('off')
    title('(b)','position',[103,-59.3],'fontweight','bold','fontSize',16);
    h=colorbar('south');
    set(h,'fontweight','bold','axislocation','out','position',[0.378616438356164,0.175845360824742,0.3,0.03]);
    set(get(h,'ylabel'), 'string','Vp (km/s)','fontsize',15,'fontweight','bold');
    hold on
    
    subplot(222)
    SZ2 = reshape(m1(:,:,Sz(n)),[Nmx,Nmy]).';
    s3  = imagesc(mx,my,SZ2);colormap(flipud(jet));
    hold on
    set(gca,'position',[0.55,0.58,0.4,0.4])
    xlabel('X(km)','fontweight','bold','fontsize',15);
    ylabel('Y(km)','fontweight','bold','fontsize',15);
    set(gca,'fontweight','bold','ydir','normal');
    colorbar('off');
    axis equal;axis([min(mx)-dmx/2 max(mx)+dmx/2 min(my)-dmy/2 max(my)+dmy/2])
    txt = {'(b)'} ;
    text(-800,5,txt,'fontsize',20) ;
    caxis([2.7,5.8]);
    
    
    subplot(224)
    SY2 = reshape(m1(:,Sy(n),:),[Nmx,Nmz]).';
    s4  = imagesc(mx,mz,SY2);colormap(flipud(jet));
    hold on
    set(gca,'position',[0.55 0.3 0.4 0.2])
    xlabel('X(km)','fontweight','bold','fontsize',15);
    ylabel('Z(km)','fontweight','bold','fontsize',15);
    set(gca,'fontweight','bold','ydir','reverse');
    colorbar('off')
    txt = {'(d)'} ;
    text(-800,5,txt,'fontsize',20);
    caxis([2.7,5.8]);
end

