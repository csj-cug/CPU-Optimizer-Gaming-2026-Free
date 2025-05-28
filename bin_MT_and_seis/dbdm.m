function dbdm1  = dbdm(omega,minv_ind,TE,ed,NiGauss,EPGauss,W,JacobDet,Edge_ind)
global mu0;
ned     = size(ed,1);
Nm      = length(minv_ind);
NGauss  = length(W);

NiGaussMat1 = zeros([size(NiGauss{1}(:,:,minv_ind)) length(NiGauss)]);

for it = 1:length(NiGauss)
    NiGaussMat1(:,:,:,it) = NiGauss{it}(:,:,minv_ind);
end

% NiGaussMat = zeros(Ne,3,NGauss,size(TE,2)); %Convert NiGauss to matrix format

NiGaussMat = permute(NiGaussMat1,[3 1 2 4]);



WRepmat = zeros(Nm,3,NGauss);

for it = 1:NGauss
    WRepmat(:,:,it)=W(it);
end

JacobDetRepmat1 = repmat(JacobDet,[1 1 3]);
JacobDetRepmat  = permute(JacobDetRepmat1,[1 3 2]);
clear JacobDetRepmat1;

ii   = zeros(12*Nm,1);
sBx  = zeros(12*Nm,1);

index = 0;
for it1= 1:12 
    ii(index+1:index+Nm) = TE(:,it1);
    jj(index+1:index+Nm) = 1:Nm;% the corresponding domain
    tmp = sum(NiGaussMat(:,:,:,it1).*JacobDetRepmat.*EPGauss.*WRepmat,3);
    tmp1 = sum(tmp,2);
    
    sBx(index+1:index+Nm) = 1i*omega*mu0*tmp1;
    
    index = index + Nm;
        
end

dbdm  = sparse(ii,jj,-sBx,ned,Nm);
dbdm1  = dbdm(Edge_ind,:);

1;
