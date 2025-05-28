function A = GravFwdTet(recpar_gra,t_xyz_center,vol)
global gam;
x = t_xyz_center(:,1);                  
y = t_xyz_center(:,2);                  
z = t_xyz_center(:,3);       
%z = flip(z);
nm = length(x);                     
comp = unique(recpar_gra(:,4));
comp = comp(comp<400);
ncomp = length(comp);
nd_all = size(recpar_gra,1);

A=zeros(nd_all,nm);

for icomp = 1:ncomp
    recpar=recpar_gra(recpar_gra(:,4)==comp(icomp),:);
    xp = recpar(:,1);                        %测点坐标
    yp = recpar(:,2);
    zp = recpar(:,3);
    nd = length(xp);
    Acomp=zeros(nd,nm);
    if(comp(icomp)==303)
        fprintf('Calculate the kernel for gz\n');
        parfor id=1:nd
            r       = sqrt((x-xp(id)).^2 + (y-yp(id)).^2 + (z-zp(id)).^2);
            Acomp(id,:) = gam*(((z-zp(id))./(r.^3)).*vol);                   %正演计算公式
        end
    elseif(comp(icomp)==331)
        fprintf('Calculate the kernel for gxx\n');
        parfor id=1:nd
            r       = sqrt((x-xp(id)).^2 + (y-yp(id)).^2 + (z-zp(id)).^2);
            Acomp(id,:) =1e4*gam*(3*(x(:)-xp(id)).^2./r.^2-1).*vol./(r.^3);
        end
    elseif(comp(icomp)==332)
        fprintf('Calculate the kernel for gxy\n');
        parfor id=1:nd
            r       = sqrt((x-xp(id)).^2 + (y-yp(id)).^2 + (z-zp(id)).^2);
            Acomp(id,:) =1e4*gam*(3*(x(:)-xp(id)).*(y(:)-yp(id))./r.^2).*vol./(r.^3);
        end
    elseif(comp(icomp)==333)
        fprintf('Calculate the kernel for gxz\n');
        parfor id=1:nd
            r       = sqrt((x-xp(id)).^2 + (y-yp(id)).^2 + (z-zp(id)).^2);
            Acomp(id,:) =1e4*gam*(3*(x(:)-xp(id)).*(z(:)-zp(id))./r.^2).*vol./(r.^3);
        end
    elseif(comp(icomp)==334)
        fprintf('Calculate the kernel for gyx\n');
        parfor id=1:nd
            r       = sqrt((x-xp(id)).^2 + (y-yp(id)).^2 + (z-zp(id)).^2);
            Acomp(id,:) =1e4*gam*(3*(x(:)-xp(id)).*(y(:)-yp(id))./r.^2).*vol./(r.^3);  
        end
    elseif(comp(icomp)==335)
        fprintf('Calculate the kernel for gyy\n');
        parfor id=1:nd
            r       = sqrt((x-xp(id)).^2 + (y-yp(id)).^2 + (z-zp(id)).^2);
            Acomp(id,:) =1e4*gam*(3*(y(:)-yp(id)).^2./r.^2-1).*vol./(r.^3);
        end
    elseif(comp(icomp)==336)
        fprintf('Calculate the kernel for gyz\n');
        parfor id=1:nd
            r       = sqrt((x-xp(id)).^2 + (y-yp(id)).^2 + (z-zp(id)).^2);
            Acomp(id,:) =1e4*gam*(3*(y(:)-yp(id)).*(z(:)-zp(id))./r.^2).*vol./(r.^3);
        end
    elseif(comp(icomp)==337)
        fprintf('Calculate the kernel for gzx\n');
        parfor id=1:nd
            r       = sqrt((x-xp(id)).^2 + (y-yp(id)).^2 + (z-zp(id)).^2);
            Acomp(id,:) =1e4*gam*(3*(x(:)-xp(id)).*(z(:)-zp(id))./r.^2).*vol./(r.^3);
        end
    elseif(comp(icomp)==338)
        fprintf('Calculate the kernel for gzy\n');
        parfor id=1:nd
            r       = sqrt((x-xp(id)).^2 + (y-yp(id)).^2 + (z-zp(id)).^2);
            Acomp(id,:) =1e4*gam*(3*(y(:)-yp(id)).*(z(:)-zp(id))./r.^2).*vol./(r.^3);
        end
    elseif(comp(icomp)==339)
        fprintf('Calculate the kernel for gzz\n');
        parfor id=1:nd
            r       = sqrt((x-xp(id)).^2 + (y-yp(id)).^2 + (z-zp(id)).^2);
            Acomp(id,:) =1e4*gam*(3*(z(:)-zp(id)).^2./r.^2-1).*vol./(r.^3);
        end
    else
        error('Invalid gravity component');
    end
    
    A(recpar_gra(:,4)==comp(icomp),:)=Acomp;
end