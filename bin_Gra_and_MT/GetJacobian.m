function [Jacob, JacobInv, JacobDet] = GetJacobian(p,t,Gauss)
Ne      = size(t,1);
NGauss  = length(Gauss)^3; %%Total Guass point in 3D;
Jacob   = cell(1,NGauss);
JacobInv= cell(1,NGauss);
JacobDet= zeros(Ne,NGauss);%The determinant of Jacobian matrix

tx      = zeros(size(t));    %% The x coordinate for each node in each element
ty      = zeros(size(t));
tz      = zeros(size(t));

chi_i = [-1 1 1 -1 -1 1 1 -1];
eta_i = [-1 -1 1 1 -1 -1 1 1];
zeta_i=[-1 -1 -1 -1 1 1 1 1];

for it = 1:size(t,2)
    tx(:,it) = p(t(:,it),1);
    ty(:,it) = p(t(:,it),2);
    tz(:,it) = p(t(:,it),3);
end

%%Define each Jacobian
for it = 1:NGauss
    Jacob{it} = zeros(3,3,Ne);
    JacobInv{it} = zeros(3,3,Ne);
end


for izeta = 1:length(Gauss)
    for ieta = 1:length(Gauss)
        for ichi = 1:length(Gauss)
            ind = (izeta-1)*length(Gauss)^2 + (ieta-1)*length(Gauss)+ichi;
            [Jacob{ind}, JacobInv{ind}, JacobDet(:,ind)] = GetJacobian_sub(tx,ty,tz,chi_i,eta_i,zeta_i,Gauss(ichi),Gauss(ieta),Gauss(izeta));
        end
    end
end


function [Jb,JbInv,JbDet] = GetJacobian_sub(tx,ty,tz,chi_i,eta_i,zeta_i,chi,eta,zeta)

nd      = size(tx,2); %number of node in each element;
Ne      = size(tx,1);
Jb      = zeros(3,3,Ne);
JbInv   = zeros(3,3,Ne); %Inverse of Jacobian
JbDet   = zeros(Ne,1);
Jb_tmp  = cell(1,nd);
for it = 1:nd
    Jb_tmp{it} = zeros(3,3,Ne);
end

for it = 1:nd
    Jb_tmp{it}(1,1,:) = chi_i(it)*(1+eta_i(it)*eta)*(1+zeta_i(it)*zeta)*tx(:,it);
    Jb_tmp{it}(1,2,:) = chi_i(it)*(1+eta_i(it)*eta)*(1+zeta_i(it)*zeta)*ty(:,it);
    Jb_tmp{it}(1,3,:) = chi_i(it)*(1+eta_i(it)*eta)*(1+zeta_i(it)*zeta)*tz(:,it);
    
    
    Jb_tmp{it}(2,1,:) = eta_i(it)*(1+chi_i(it)*chi)*(1+zeta_i(it)*zeta)*tx(:,it);
    Jb_tmp{it}(2,2,:) = eta_i(it)*(1+chi_i(it)*chi)*(1+zeta_i(it)*zeta)*ty(:,it);
    Jb_tmp{it}(2,3,:) = eta_i(it)*(1+chi_i(it)*chi)*(1+zeta_i(it)*zeta)*tz(:,it);
    
    
    Jb_tmp{it}(3,1,:) = zeta_i(it)*(1+chi_i(it)*chi)*(1+eta_i(it)*eta)*tx(:,it);
    Jb_tmp{it}(3,2,:) = zeta_i(it)*(1+chi_i(it)*chi)*(1+eta_i(it)*eta)*ty(:,it);
    Jb_tmp{it}(3,3,:) = zeta_i(it)*(1+chi_i(it)*chi)*(1+eta_i(it)*eta)*tz(:,it);
end

for it = 1:nd
    Jb = Jb+Jb_tmp{it}/8;
end

parfor it = 1:Ne
    JbInv(:,:,it) = inv(Jb(:,:,it));
end

JbDet = squeeze(Jb(1,1,:).*Jb(2,2,:).*Jb(3,3,:) - Jb(1,1,:).*Jb(2,3,:).*Jb(3,2,:)...
       +Jb(1,2,:).*Jb(2,3,:).*Jb(3,1,:) - Jb(1,2,:).*Jb(2,1,:).*Jb(3,3,:)...
       +Jb(1,3,:).*Jb(2,1,:).*Jb(3,2,:) - Jb(1,3,:).*Jb(2,2,:).*Jb(3,1,:));
    
% JbDet1=zeros(size(JbDet));
% for it = 1:Ne
%     JbDet1(it)=det(Jb(:,:,it));
% end
% dif = JbDet-JbDet1;
% unique(dif)
% 
% 1;







