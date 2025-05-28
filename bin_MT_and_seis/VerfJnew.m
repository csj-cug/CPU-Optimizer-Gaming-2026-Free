function VerfJnew(J,dpr,freq,ns,recpar,par,omega,np,npo,MeshTop,bnd,hl,sl,al,sig0x,sig0y,sig0z,sigx,sigy,...
           sigz,Jacob,JacobInv,JacobDet,Gauss3DW,Gauss,EPGauss,compress,symFlag,fileformat,Num_Comm,Q_Ex,Q_Ey,Q_Ez,Q_Hx,Q_Hy,Q_Hz,mpirun,xa,ya,za,minv_ind,Xtop,Ytop,Ztop,inv_flag) 
%%Verify sensistivity matrix
%Using perturbance method

t_xyz_center = MeshTop.t_xyz_center;

%find anamous body doamin
anamous_ind0  = find(t_xyz_center(:,1)>=xa{1}(1) & t_xyz_center(:,1)<=xa{1}(2) & t_xyz_center(:,2)>=ya{1}(1)...
                 & t_xyz_center(:,2)<=ya{1}(2) & t_xyz_center(:,3)>=za{1}(1) & t_xyz_center(:,3)<=za{1}(2));   
% anamous_ind   = anamous_ind0(100);             
sigx(anamous_ind0) = sigx(anamous_ind0).*(1+0.01);
sigy(anamous_ind0) = sigy(anamous_ind0).*(1+0.01);
sigz(anamous_ind0) = sigz(anamous_ind0).*(1+0.01);

d  = FEM_forward(freq,ns,recpar,par,omega,np,npo,MeshTop,bnd,hl,sl,al,sig0x,sig0y,sig0z,sigx,sigy,...
       sigz,Jacob,JacobInv,JacobDet,Gauss3DW,Gauss,EPGauss,compress,symFlag,fileformat,Num_Comm,Q_Ex,Q_Ey,Q_Ez,Q_Hx,Q_Hy,Q_Hz,mpirun,Xtop,Ytop,Ztop,inv_flag);

Jpert = (d - dpr)./0.01;

%Sensistivity matrix from reciprocity method
Jrecp = zeros(size(J,1),1);
for it = 1:length(anamous_ind0)
    Jrecp = Jrecp + J(:,minv_ind == anamous_ind0(it));
end

% %
figure
subplot(2,1,1)
plot(real(Jrecp),'ro-');
hold on
plot(real(Jpert),'b*-');
hold off
title('Real part of sensisvity')
xlabel('Observation data');
ylabel('Amplitude');
legend('Adjoint method','Perturbance method');

subplot(2,1,2)
plot(imag(Jrecp),'ro-');
hold on
plot(imag(Jpert),'b*-');
hold off
xlabel('Observation data');
ylabel('Amplitude');
title('Imag part of sensisvity')
legend('Adjoint method','Perturbance method');
