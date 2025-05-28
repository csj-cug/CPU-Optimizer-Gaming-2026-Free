function b = getB_edge_New(omega,sig0x,sig0y,sig0z,sigx,sigy,sigz,TE,ed,NiGauss,EPGauss,W,JacobDet)
global mu0;
ned     = size(ed,1);
Ne      = size(TE,1);
NGauss  = length(W);

sigax     = sigx-sig0x;%%Anomalous conductivity
sigay     = sigy-sig0y;
sigaz     = sigz-sig0z;

siga     = zeros(Ne,3,NGauss); % 
siga(:,1,:) = repmat(sigax,[1,1,NGauss]);
siga(:,2,:) = repmat(sigay,[1,1,NGauss]);
siga(:,3,:) = repmat(sigaz,[1,1,NGauss]);


NiGaussMat1 = zeros([size(NiGauss{1}) length(NiGauss)]);

for it = 1:length(NiGauss)
    NiGaussMat1(:,:,:,it) = NiGauss{it};
end

% NiGaussMat = zeros(Ne,3,NGauss,size(TE,2)); %Convert NiGauss to matrix format

NiGaussMat = permute(NiGaussMat1,[3 1 2 4]);



WRepmat = zeros(Ne,3,NGauss);

for it = 1:NGauss
    WRepmat(:,:,it)=W(it);
end


JacobDetRepmat1 = repmat(JacobDet,[1 1 3]);
JacobDetRepmat  = permute(JacobDetRepmat1,[1 3 2]);
clear JacobDetRepmat1;


ii   = zeros(12*Ne,1);
sBx  = zeros(12*Ne,1);

index = 0;
for it1= 1:12 
    ii(index+1:index+Ne) = TE(:,it1);
    tmp = sum(NiGaussMat(:,:,:,it1).*JacobDetRepmat.*siga.*EPGauss.*WRepmat,3);
    tmp1 = sum(tmp,2);
    
    sBx(index+1:index+Ne) = 1i*omega*mu0*tmp1;
    
    index = index + Ne;
        
end


b  = sparse(ii,1,-sBx,ned,1);

b = full(b);

1;
