function EPGauss = GetXYZGauss(p,t,L,epxnode,epynode,epznode)

Ne      = size(t,1);
NGauss  = size(L,1);

Xnode   = [p(t(:,1),1), p(t(:,2),1), p(t(:,3),1),p(t(:,4),1),p(t(:,5),1), p(t(:,6),1), p(t(:,7),1),p(t(:,8),1)];
Ynode   = [p(t(:,1),2), p(t(:,2),2), p(t(:,3),2),p(t(:,4),2),p(t(:,5),2), p(t(:,6),2), p(t(:,7),2),p(t(:,8),2)];
Znode   = [p(t(:,1),3), p(t(:,2),3), p(t(:,3),3),p(t(:,4),3),p(t(:,5),3), p(t(:,6),3), p(t(:,7),3),p(t(:,8),3)];




TE_epx = zeros(size(t));
TE_epy = zeros(size(t));
TE_epz = zeros(size(t));

for it = 1:8
    TE_epx(:,it)=epxnode(t(:,it));
    TE_epy(:,it)=epynode(t(:,it));
    TE_epz(:,it)=epznode(t(:,it));
end

EPGauss       = zeros(Ne,3,NGauss);

for it=1:NGauss
    LRepmat = repmat(L(it,:),Ne,1);
    
    EPGauss(:,1,it) = sum(TE_epx.*LRepmat,2);
    EPGauss(:,2,it) = sum(TE_epy.*LRepmat,2);
    EPGauss(:,3,it) = sum(TE_epz.*LRepmat,2);
end
