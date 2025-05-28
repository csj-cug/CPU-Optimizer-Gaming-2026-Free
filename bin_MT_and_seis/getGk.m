function Gk1=getGk(omega1,omega,Dm3D,TE,ed,M,Edge_ind,minv_ind)
ned     = size(ed,1);
Nm      = size(TE,1);

ii       = Dm3D(minv_ind,1:12*12);
jj       = Dm3D(minv_ind,12*12 + 1:2*12*12);
sF_Gauss = Dm3D(minv_ind,2*12*12 + 1:3*12*12);

ii1       = zeros(12*12*Nm,1);
jj1       = zeros(12*12*Nm,1);
kk1       = zeros(12*12*Nm,1);
sF_Gauss2 = zeros(12*12*Nm,1);
index = 0;
for it = 1:12
    for jt = 1:12
        ii1(index+1:index+Nm,1) = ii(:,(it-1)*12 + jt);
        jj1(index+1:index+Nm,1) = jj(:,(it-1)*12 + jt);
        kk1(index+1:index+Nm,1) = 1:Nm;
        sF_Gauss1 = sF_Gauss(:,(it-1)*12 + jt);
        sF_Gauss2(index+1:index+Nm,1) = omega/omega1*sF_Gauss1.*M(ii1(index+1:index+Nm,1));
        index = index + Nm;
    end
end
Gk = sparse(jj1,kk1,sF_Gauss2,ned,Nm);
Gk1=Gk(Edge_ind,:);
