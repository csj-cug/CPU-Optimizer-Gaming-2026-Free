function [K, K1, Ni, CurlNi] = getK_Brick_Edge_Jacobian(omega,sigx,sigy,sigz,TE,TE_length,ed,Jacob,JacobInv,JacobDet,Gauss,Gauss3DW)
global mu0;
ned     = size(ed,1);
Ne      = size(TE,1);
nGauss  = length(Gauss);
%% Define the center of each edge in \xi-\eta-\zeta domain
chi_edge  = [0 0 0 0 -1 -1 1 1 -1 1 -1 1];
eta_edge  = [-1 1 -1 1 0 0 0 0 -1 -1 1 1];
zeta_edge = [-1 -1 1 1 -1 1 -1 1 0 0 0 0];

%%The gradient defined on each Gauss point
Grad_chi     = zeros(3,nGauss^3,Ne);
Grad_eta     = zeros(3,nGauss^3,Ne);
Grad_zeta    = zeros(3,nGauss^3,Ne);

for it = 1:nGauss^3
    Grad_chi(:,it,:)  = JacobInv{it}(:,1,:);
    Grad_eta(:,it,:)  = JacobInv{it}(:,2,:);
    Grad_zeta(:,it,:) = JacobInv{it}(:,3,:);  
end

%%Define the edge-based shape function it its curl
Ni     = cell(1,12);
CurlNi = cell(1,12);

for it = 1:12
    Ni{it}     = zeros(3,nGauss^3,Ne);
    CurlNi{it} = zeros(3,nGauss^3,Ne);
end


chi_gauss   = zeros(1,nGauss^3);
eta_gauss   = zeros(1,nGauss^3);
zeta_gauss  = zeros(1,nGauss^3);
for izeta  = 1:length(Gauss)
    for ieta = 1:length(Gauss)
        for ichi = 1:length(Gauss)
            ind = (izeta-1)*length(Gauss)^2 + (ieta-1)*length(Gauss)+ichi;
            chi_gauss(ind)  = Gauss(ichi);
            eta_gauss(ind)  = Gauss(ieta);
            zeta_gauss(ind) = Gauss(izeta);
            
            
        end
    end
end



for it =1:12
    Lg = zeros(3,nGauss^3,Ne);
    for ii = 1:3
        for jj = 1:nGauss^3
            Lg(ii,jj,:) = TE_length(:,it);
        end
    end
    tmp = [];
    if it>=1 && it<=4
        tmp = Lg.*(eta_edge(it)*Grad_eta + zeta_edge(it)*Grad_zeta + ...
            eta_edge(it)*zeta_edge(it)*repmat(eta_gauss,[3,1,Ne]).*Grad_zeta + eta_edge(it)*zeta_edge(it)*repmat(zeta_gauss,[3,1,Ne]).*Grad_eta)/8; 
        CurlNi{it} = cross(tmp,Grad_chi,1);
        Ni{it}     = Lg.*((1+eta_edge(it)*repmat(eta_gauss,[3,1,Ne])).*(1+zeta_edge(it)*repmat(zeta_gauss,[3,1,Ne])).*Grad_chi)/8;
             
    elseif it>=5 && it <=8
        tmp = Lg.*(chi_edge(it)*Grad_chi + zeta_edge(it)*Grad_zeta + ...
            chi_edge(it)*zeta_edge(it)*repmat(chi_gauss,[3,1,Ne]).*Grad_zeta + chi_edge(it)*zeta_edge(it)*repmat(zeta_gauss,[3,1,Ne]).*Grad_chi)/8; 
        CurlNi{it} = cross(tmp,Grad_eta,1);
        Ni{it}     = Lg.*((1+chi_edge(it)*repmat(chi_gauss,[3,1,Ne])).*(1+zeta_edge(it)*repmat(zeta_gauss,[3,1,Ne])).*Grad_eta)/8;
        
    elseif it>=9 && it<=12
        tmp = Lg.*(chi_edge(it)*Grad_chi + eta_edge(it)*Grad_eta + ...
            chi_edge(it)*eta_edge(it)*repmat(chi_gauss,[3,1,Ne]).*Grad_eta + chi_edge(it)*eta_edge(it)*repmat(eta_gauss,[3,1,Ne]).*Grad_chi)/8; 
        CurlNi{it} = cross(tmp,Grad_zeta,1);
        Ni{it}     = Lg.*((1+chi_edge(it)*repmat(chi_gauss,[3,1,Ne])).*(1+eta_edge(it)*repmat(eta_gauss,[3,1,Ne])).*Grad_zeta)/8;
    end
end


sig = zeros(3,nGauss^3,Ne);

for it =1:nGauss^3
    sig(1,it,:)=sigx;
    sig(2,it,:)=sigy;
    sig(3,it,:)=sigz;
end

sigc = zeros(3,nGauss^3,Ne);     %to compute the G,DA/Dm
for it =1:nGauss^3
    sigc(1,it,:)=1;
    sigc(2,it,:)=1;
    sigc(3,it,:)=1;
end



ned     = size(ed,1);
Ne      = size(TE,1);
ii = zeros(12*12*Ne,1);
jj = zeros(12*12*Ne,1);
sE_Gauss  = zeros(12*12*Ne,1);
sF_Gauss  = zeros(12*12*Ne,1);
sF_Gauss1 = zeros(12*12*Ne,1);
index = 0;




for it1= 1:12
    for jt1 = 1:12
        ii(index+1:index+Ne) = TE(:,it1);
        jj(index+1:index+Ne) = TE(:,jt1);
        
%         add the weighting inside, now assume weight is 1;
        
        CurlNi_Dot = squeeze(dot(CurlNi{it1},CurlNi{jt1},1))';
        Ni_Dot     = 1i*squeeze(dot(Ni{it1},(omega*mu0*sig).*Ni{jt1}))';
        Ni_Dot1    = 1i*squeeze(dot(Ni{it1},(omega*mu0*sigc).*Ni{jt1}))';
        
        sE_Gauss(index+1:index+Ne) = sum(CurlNi_Dot.*JacobDet.*repmat(Gauss3DW,Ne,1),2);
        sF_Gauss(index+1:index+Ne) = sum(Ni_Dot.*JacobDet.*repmat(Gauss3DW,Ne,1),2);
        sF_Gauss1(index+1:index+Ne) = sum(Ni_Dot1.*JacobDet.*repmat(Gauss3DW,Ne,1),2);
        index = index + Ne;                       
    end
end
K = sparse(ii,jj,sE_Gauss+sF_Gauss,ned,ned);

%rewrite ii, jj, sF_Gauss1 each element
ii1        = reshape(ii,Ne,12*12);
jj1        = reshape(jj,Ne,12*12);
sF_Gauss2 = reshape(sF_Gauss1,Ne,12*12);
K1 = [ii1,jj1,sF_Gauss2];

