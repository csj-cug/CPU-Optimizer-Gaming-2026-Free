function [sig0x, sig0y, sig0z, sigx, sigy, sigz] = SetSig_New(t_xyz_center,t_xyz_center_reg,xa,ya,za,sigax,sigay,sigaz,hl,sl,al,loadSig,AnoGridReg,viewflag)
global sig_air;
Ne      = size(t_xyz_center,1);
slv     = sl./(al.^2);       %Vector of vertical conductivity
sig0x   = zeros(Ne,1);       %Background conductivity
sig0y   = zeros(Ne,1);
sig0z   = zeros(Ne,1);
sigx    = zeros(Ne,1);       %Total conductivity
sigy    = zeros(Ne,1);
sigz    = zeros(Ne,1);
nlayer  = length(sl);        %Number of layers
na      = length(xa);  %Number of anomalous domain

x  = unique(t_xyz_center_reg(:,1));
y  = unique(t_xyz_center_reg(:,2));
nx = length(x);
ny = length(y);
nz = Ne/(nx*ny);

Xc_reg=reshape(t_xyz_center_reg(:,1),nx,ny,nz);
Yc_reg=reshape(t_xyz_center_reg(:,2),nx,ny,nz);
Zc_reg=reshape(t_xyz_center_reg(:,3),nx,ny,nz);

if loadSig==0
    if na~=length(sigax)
        error('Number of anomlous domain does not match with the sigax');
    end
end

sig0x((t_xyz_center(:,3)<0)) = sig_air;
sig0y((t_xyz_center(:,3)<0)) = sig_air;
sig0z((t_xyz_center(:,3)<0)) = sig_air;


sigx((t_xyz_center_reg(:,3)<0)) = sig_air;
sigy((t_xyz_center_reg(:,3)<0)) = sig_air;
sigz((t_xyz_center_reg(:,3)<0)) = sig_air;

for ilayer = 1:nlayer-1
    if ilayer == 1
        z1 = 0;
        z2 = hl(1);
    else
        z1 = sum(hl(1:ilayer-1));
        z2 = sum(hl(1:ilayer));
    end
    sig0x((t_xyz_center(:,3)>z1)&(t_xyz_center(:,3)<z2)) = sl(ilayer); %%Set the background conductivity for each cell according to its true coordinate for the cell center
    sig0y((t_xyz_center(:,3)>z1)&(t_xyz_center(:,3)<z2)) = sl(ilayer);
    sig0z((t_xyz_center(:,3)>z1)&(t_xyz_center(:,3)<z2)) = slv(ilayer);
    
    sigx((t_xyz_center_reg(:,3)>z1)&(t_xyz_center_reg(:,3)<z2)) = sl(ilayer); %%Set the true conductivity for each cell according to its reference cell center (regular grid)
    sigy((t_xyz_center_reg(:,3)>z1)&(t_xyz_center_reg(:,3)<z2)) = sl(ilayer);
    sigz((t_xyz_center_reg(:,3)>z1)&(t_xyz_center_reg(:,3)<z2)) = slv(ilayer);
end

z3    = sum(hl);
sig0x(t_xyz_center(:,3)>z3) = sl(end); %%Set the background conductivity for the last layer
sig0y(t_xyz_center(:,3)>z3) = sl(end);
sig0z(t_xyz_center(:,3)>z3) = slv(end);

sigx(t_xyz_center_reg(:,3)>z3) = sl(end);%%Set the true conductivity for the last layer
sigy(t_xyz_center_reg(:,3)>z3) = sl(end);
sigz(t_xyz_center_reg(:,3)>z3) = slv(end);

%Check if sig0x = sig0y
dif0 = unique(sig0x-sig0y);
dif = unique(sigx-sigy);
if length(dif0)~=1 || length(dif)~=1
    error('Horizontal background conductivity is not equal to vertical one');
else
    if dif0~=0 || dif~=0
        error('Horizontal background conductivity is not equal to vertical one');
    end
end


Sigx = reshape(sigx,nx,ny,nz); %covert sigx to 3D matrix
Sigy = reshape(sigy,nx,ny,nz);
Sigz = reshape(sigz,nx,ny,nz);


if loadSig
    load sigbody;
end

for ia = 1:na
    if AnoGridReg==0 %% define the anomalous domain based on the true grid after distortion
        if loadSig == 0 %define the anomalous conductivity in the Input.m file
            ind = find(t_xyz_center(:,1)>=xa{ia}(1)&t_xyz_center(:,1)<=xa{ia}(2)&...
                t_xyz_center(:,2)>=ya{ia}(1)&t_xyz_center(:,2)<=ya{ia}(2)&...
                t_xyz_center(:,3)>=za{ia}(1)&t_xyz_center(:,3)<=za{ia}(2));
            sigx(ind)=sigax(ia);
            sigy(ind)=sigay(ia);
            sigz(ind)=sigaz(ia);
        else
            error('Define anomalous domain based on distorted grid, the loadSig option is not working now; Please define the anomalous conductivity in the Input.m file');
        end
        
    elseif AnoGridReg==1 %% define the anomalous domain based on the regular grid before distortion
        if loadSig == 0 %define the anomalous conductivity in the Input.m file
            %             ind = find(t_xyz_center_reg(:,1)>=xa{ia}(1)&t_xyz_center_reg(:,1)<=xa{ia}(2)&...
            %                 t_xyz_center_reg(:,2)>=ya{ia}(1)&t_xyz_center_reg(:,2)<=ya{ia}(2)&...
            %                 t_xyz_center_reg(:,3)>=za{ia}(1)&t_xyz_center_reg(:,3)<=za{ia}(2));
            %
            %             sigx(ind)=sigax(ia);
            %             sigy(ind)=sigay(ia);
            %             sigz(ind)=sigaz(ia);
            
            
            
            indx = find(unique(Xc_reg)>=xa{ia}(1)&unique(Xc_reg)<=xa{ia}(2));
            indy = find(unique(Yc_reg)>=ya{ia}(1)&unique(Yc_reg)<=ya{ia}(2));
            indz = find(unique(Zc_reg)>=za{ia}(1)&unique(Zc_reg)<=za{ia}(2));
            nx_ano = length(indx); ny_ano = length(indy); nz_ano = length(indz);
            Ne_ano = nx_ano*ny_ano*nz_ano; %number of cells in the anomalous domain
            
            Sigx(indx,indy,indz) = sigax(ia);
            Sigy(indx,indy,indz) = sigay(ia);
            Sigz(indx,indy,indz) = sigaz(ia);
            
        else %load the conductivity for the anomaluos domain
            if size(sigbody{ia},2) ~= 3
                error('The number of columns for sigbody should be 3 times numer of anomalous domain');
            end
            indx = find(unique(Xc_reg)>=xa{ia}(1)&unique(Xc_reg)<=xa{ia}(2));
            indy = find(unique(Yc_reg)>=ya{ia}(1)&unique(Yc_reg)<=ya{ia}(2));
            indz = find(unique(Zc_reg)>=za{ia}(1)&unique(Zc_reg)<=za{ia}(2));
            nx_ano = length(indx); ny_ano = length(indy); nz_ano = length(indz);
            Ne_ano = nx_ano*ny_ano*nz_ano; %number of cells in the anomalous domain
            
            if size(sigbody{ia},1) ==1
                Sigx(indx,indy,indz) = sigbody{ia}(1);
                Sigy(indx,indy,indz) = sigbody{ia}(2);
                Sigz(indx,indy,indz) = sigbody{ia}(3);
            elseif size(sigbody{ia},1) == Ne_ano
                Sigx(indx,indy,indz) = reshape(sigbody{ia}(:,1),nx_ano,ny_ano,nz_ano);
                Sigy(indx,indy,indz) = reshape(sigbody{ia}(:,2),nx_ano,ny_ano,nz_ano);
                Sigz(indx,indy,indz) = reshape(sigbody{ia}(:,3),nx_ano,ny_ano,nz_ano);
            else
                error('size of sigbody does not match the number of cells in the anomalous domain');
            end
            
        end
        
        sigx = Sigx(:);
        sigy = Sigy(:);
        sigz = Sigz(:);
        
    else
        error('AnoGridReg should be either 0 or 1');
    end
end

if viewflag
    show_model_hex(t_xyz_center,sig0x,0,[],'\sigma_0_x',1);
    show_model_hex(t_xyz_center,sig0y,0,[],'\sigma_0_y',1);
    show_model_hex(t_xyz_center,sig0z,0,[],'\sigma_0_z',1);
    show_model_hex(t_xyz_center,sigx,0,[],'\sigma_x',1);
    show_model_hex(t_xyz_center,sigy,0,[],'\sigma_y',1);
    show_model_hex(t_xyz_center,sigz,0,[],'\sigma_z',1);
    show_model_hex(t_xyz_center,sigx-sig0x,0,[],'\sigma_a_h',1); cx=caxis;caxis([-max(abs(cx)) max(abs(cx))]);
    show_model_hex(t_xyz_center,sigz-sig0z,0,[],'\sigma_a_v',1); cx=caxis;caxis([-max(abs(cx)) max(abs(cx))]);
end
% 
% 
% if viewflag
%     xi = linspace(-1500,1500,40);
%     yi = linspace(-1500,1500,40);
%     zi = linspace(-500,4000,100);
%     xi = [-1000:100:-500 -480:20:480 500:100:1000];
%     yi = xi;
%     zi = [-500:100:-100 -50 -20 20 50 100 200:100:800 850 910:10:1100 1150:50:1500 1600:100:2000 2020:20:2100 2120 2200:100:2500];
%     [xi, yi, zi] = meshgrid(xi,yi,zi);
% 
%     vi0x = griddata(t_xyz_center(:,1),t_xyz_center(:,2),t_xyz_center(:,3),sig0x,xi,yi,zi,'nearest');
%     vi0y = griddata(t_xyz_center(:,1),t_xyz_center(:,2),t_xyz_center(:,3),sig0y,xi,yi,zi,'nearest');
%     vi0z = griddata(t_xyz_center(:,1),t_xyz_center(:,2),t_xyz_center(:,3),sig0z,xi,yi,zi,'nearest');
% 
%     vix  = griddata(t_xyz_center(:,1),t_xyz_center(:,2),t_xyz_center(:,3),sigx,xi,yi,zi,'nearest');
%     viy  = griddata(t_xyz_center(:,1),t_xyz_center(:,2),t_xyz_center(:,3),sigy,xi,yi,zi,'nearest');
%     viz  = griddata(t_xyz_center(:,1),t_xyz_center(:,2),t_xyz_center(:,3),sigz,xi,yi,zi,'nearest');
% 
%     figure;
%     slice(xi,yi,zi,vi0x,[0],[0],[]);
%     set(gca,'zdir','reverse');
%     colorbar;
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     shading flat;
%     alpha(.8);
%     title('\sigma_0_x');
% 
% 
% 
%     figure;
%     slice(xi,yi,zi,vi0y,[0],[0],[]);
%     set(gca,'zdir','reverse');
%     colorbar;
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     shading flat;
%     alpha(.8);
%     title('\sigma_0_y');
% 
% 
%     figure;
%     slice(xi,yi,zi,vi0z,[0],[0],[]);
%     set(gca,'zdir','reverse');
%     colorbar;
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     shading flat;
%     alpha(.8);
%     title('\sigma_0_z');
% 
% 
%     figure;
%     slice(xi,yi,zi,vix,[0],[0],[]);
%     set(gca,'zdir','reverse');
%     colorbar;
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     shading flat;
%     alpha(.8);
%     title('\sigma_x');
% 
% 
% 
%     figure;
%     slice(xi,yi,zi,viy,[0],[0],[]);
%     set(gca,'zdir','reverse');
%     colorbar;
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     shading flat;
%     alpha(.8);
%     title('\sigma_y');
% 
% 
%     figure;
%     slice(xi,yi,zi,viz,[0],[0],[]);
%     set(gca,'zdir','reverse');
%     colorbar;
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     shading flat;
%     alpha(.8);
%     title('\sigma_z');
% end

