function show_model_hex(t_xyz_center,sig,xslice,yslice,fig_title,flipcolormap)
%%show the x,y section of the hexehedral model
%t_xyz_center: cell center for the hex mesh
%sig: conductivity vector for each cell,
%xslice: x coordinate of the x xslice
%zslice: y coordinate of the y xslice
%fig_title: a string for part of the figure title
%flipcolormap: 0: use jet, 1 or other value:flipud(jet);
Ne = size(t_xyz_center,1);
x  = unique(t_xyz_center(:,1));
y  = unique(t_xyz_center(:,2));
nx = length(x);
ny = length(y);
nz = Ne/(nx*ny);




Xc=reshape(t_xyz_center(:,1),nx,ny,nz);
Yc=reshape(t_xyz_center(:,2),nx,ny,nz);
Zc=reshape(t_xyz_center(:,3),nx,ny,nz);
Sig=reshape(sig,nx,ny,nz);

% if ~isempty(yslice)
%     for iy=1:length(yslice)
%         ind = find(abs(yslice(iy)-y)==min(abs(yslice(iy)-y)));
%         ind = ind(1);
%         figure;
%         hp=pcolor(squeeze(Xc(:,ind,:)),squeeze(Zc(:,ind,:)),squeeze(Sig(:,ind,:)));
%         set(gca,'ydir','reverse');
%         h=colorbar;
%         if flipcolormap
%             colormap(color); 
%         end
%         set(hp,'EdgeAlpha',.25);
%         set(h,'fontweight','bold');
%         set(gca,'fontweight','bold');
%         xlabel('x(m)');
%         ylabel('z(m)');
%         axis equal;axis([min(x) max(x) min(t_xyz_center(:,3)) max(t_xyz_center(:,3))]);
%         set(get(h,'xlabel'),'string','S/m','fontweight','bold');
%         title([fig_title ', y=' num2str(yslice(iy))]);
%     end
% end




if ~isempty(yslice)
    for iy=1:length(yslice)
        ind = find(abs(yslice(iy)-y)==min(abs(yslice(iy)-y)));
        ind = ind(1);
        figure;
        hp=pcolor(squeeze(Xc(:,ind,:)),squeeze(Zc(:,ind,:)),squeeze(Sig(:,ind,:)));
        set(gca,'ydir','reverse');
        h=colorbar;
        if flipcolormap
            colormap(flipud(jet)); 
        end
        set(hp,'EdgeAlpha',.25);
        set(h,'fontweight','bold');
        set(gca,'fontweight','bold');
        xlabel('x(m)');
        ylabel('z(m)');
        axis equal;axis([min(x) max(x) min(t_xyz_center(:,3)) max(t_xyz_center(:,3))]);
        set(get(h,'xlabel'),'string','S/m','fontweight','bold');
        title([fig_title ', y=' num2str(yslice(iy))]);
    end
end





if ~isempty(xslice)
    for ix=1:length(xslice)
        ind = find(abs(xslice(ix)-x)==min(abs(xslice(ix)-x)));
        ind = ind(1);
        figure;
        hp=pcolor(squeeze(Yc(ind,:,:)),squeeze(Zc(ind,:,:)),squeeze(Sig(ind,:,:)));
        set(gca,'ydir','reverse');
        h=colorbar;
        if flipcolormap
            colormap(flipud(jet)); 
        end
        set(hp,'EdgeAlpha',.25);
        set(h,'fontweight','bold');
        set(gca,'fontweight','bold');
        xlabel('y(m)');
        ylabel('z(m)');
        axis equal;axis([min(y) max(y) min(t_xyz_center(:,3)) max(t_xyz_center(:,3))]);
        set(get(h,'xlabel'),'string','S/m','fontweight','bold');
        title([fig_title ', x=' num2str(xslice(ix))]);
    end
end
