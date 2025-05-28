function [d,JZ] =FEM_forwardJ(nrec,freq,ns,omega,recpar,par,np,npo,MeshTop,bnd,hl,sl,al,sig0x,sig0y,sig0z,sigx...
        ,sigy,sigz,Jacob,JacobInv,JacobDet,Gauss3DW,Gauss,EPGauss,compress,symFlag,fileformat,Num_Comm,Q_Ex,Q_Ey,Q_Ez,Q_Hx,Q_Hy,Q_Hz,mpirun,Edge_ind,minv_ind,Xtop,Ytop,Ztop,inv_flag)

global mu0;

nfreq=length(freq);
[K, Dm3D, Ni,~] = getK_Brick_Edge_Jacobian(omega(1),sigx,sigy,sigz,...
    MeshTop.TE,MeshTop.TE_length,MeshTop.ed,Jacob,JacobInv,JacobDet,Gauss,Gauss3DW);
b = cell(nfreq,ns);

for ifreq = 1:nfreq
    for isrc = 1:ns
        b{ifreq,isrc} = getB_edge_New(omega(ifreq),sig0x,sig0y,sig0z,sigx...
            ,sigy,sigz,MeshTop.TE,MeshTop.ed,Ni,EPGauss{ifreq,isrc},Gauss3DW,JacobDet);
    end
end
Bdirc          = getbnd(bnd);                                              % Get boundary condition matrix
[Kglob,bglob]  = DiricBndNew_MultSrc_MultFreq(K,b,Bdirc,freq,nfreq,ns);    % Apply Dirichlet boundary condition
[~, nfreq_per_comm] = PardisoInputGeneralMultRHSMultFreq_Split(Kglob,bglob,1,symFlag,fileformat,compress,freq,Num_Comm);

mpibin = '/public/home/caihongzhu/zhuhao/MT3D/MT3D_structured_Inv/pardiso_scratch_MpiSplit/intel_intelmpi_lp64_intel64_so';
%/public/home/caihongzhu/zhuhao/Joint/Joint_structured_Inv/pardiso_scratch_MpiSplit/intel_intelmpi_lp64_intel64_so


    if ispc
        dirsep  = '\';
    else
	dirsep  = '/';
    end
    MPI_    = mpirun;
    MPIBIN_ = mpibin;
    MPIOPS_ = ['export OMP_NUM_THREADS=' num2str(npo) ' && ' MPI_  ' -np ' num2str(np) ' -machinefile nd' ' ']; 
    SOLVEROUT = 'pardiso.out';
    SOLVERSTR = 'MKL_Pardiso';
    SOLVER = 'pardiso_complex_driver.exe'; 

if strcmp(MPIBIN_,'')
    MPIDIR_ = '';
else
    MPIDIR_ = [MPIBIN_ dirsep];
end

CMD_ = [MPIOPS_ MPIDIR_ SOLVER ' >& ' SOLVEROUT];
disp('Calling external direct solver, this may take a while ...');
disp(CMD_);
system(CMD_);


M_Comm = cell(Num_Comm,1);        
for icomm = 1:Num_Comm
    M_Comm{icomm} = ReadSol_MultRHS_MultFreq(ns*nfreq_per_comm(icomm),['Output_Comm_' num2str(icomm-1) '.h5'],fileformat,compress);
end
M = [];
for icomm = 1:Num_Comm
    M = [M,M_Comm{icomm}];
end

%Get electromagnetic fields
intoutFEM = [];

for ifreq = 1:nfreq
    for isrc = 1:ns
        intoutFEMtmp = getfield_new(M(:,(ifreq-1)*ns+isrc),Q_Ex,Q_Ey,Q_Ez,Q_Hx,Q_Hy,Q_Hz,recpar,hl,sl,al,par{isrc},freq(ifreq),omega(ifreq),isrc,Xtop,Ytop,Ztop,MeshTop.p);
        intoutFEM = [intoutFEM;intoutFEMtmp];
    end
end

%Get impedances
[inptdata,~] = GetImpedanceFull(intoutFEM);

switch inv_flag
    case 1
        inptZ = [];
        for ifreq = 1:nfreq
            inptf   = inptdata(inptdata(:,5) == freq(ifreq),:);
            inptZxy = inptf(inptf(:,4) == 1,:);
            inptZyx = inptf(inptf(:,4) == 2,:);
            inptZ   = [inptZ;inptZxy;inptZyx];
        end
        npara = 4;
    case 2
        inptZ = [];
        for ifreq = 1:nfreq
            inptf   = inptdata(inptdata(:,5) == freq(ifreq),:);
            inptZxy = inptf(inptf(:,4) == 1,:);
            inptZyx = inptf(inptf(:,4) == 2,:);
            inptZxx = inptf(inptf(:,4) == 3,:);
            inptZyy = inptf(inptf(:,4) == 4,:);
            inptZ   = [inptZ;inptZxy;inptZyx;inptZxx;inptZyy];
        end
        npara = 4;
    case 3
        inptZ = [];
        for ifreq = 1:nfreq
            inptf   = inptdata(inptdata(:,5) == freq(ifreq),:);
            inptTx  = inptf(inptf(:,4) == 5,:);
            inptTy  = inptf(inptf(:,4) == 6,:);
            inptZ   = [inptZ;inptTx;inptTy];
        end   
        npara = 3;
    case 4
        inptZ = [];
        for ifreq = 1:nfreq
            inptf   = inptdata(inptdata(:,5) == freq(ifreq),:);
            inptZxy = inptf(inptf(:,4) == 1,:);
            inptZyx = inptf(inptf(:,4) == 2,:);
            inptTx  = inptf(inptf(:,4) == 5,:);
            inptTy  = inptf(inptf(:,4) == 6,:);
            inptZ   = [inptZ;inptZxy;inptZyx;inptTx;inptTy];
        end
        npara = 5;
    case 5
        inptZ = [];
        for ifreq = 1:nfreq
            inptf   = inptdata(inptdata(:,5) == freq(ifreq),:);
            inptZxy = inptf(inptf(:,4) == 1,:);
            inptZyx = inptf(inptf(:,4) == 2,:);
            inptZxx = inptf(inptf(:,4) == 3,:);
            inptZyy = inptf(inptf(:,4) == 4,:);
            inptTx  = inptf(inptf(:,4) == 5,:);
            inptTy  = inptf(inptf(:,4) == 6,:);
            inptZ   = [inptZ;inptZxy;inptZyx;inptZxx;inptZyy;inptTx;inptTy];
        end
        npara = 5;
    case 6
        inptZ = [];
        for ifreq = 1:nfreq
            inptf    = inptdata(inptdata(:,5) == freq(ifreq),:);
            inptPTxx = inptf(inptf(:,4) == 7,:);
            inptPTyy = inptf(inptf(:,4) == 8,:);
            inptPTxy = inptf(inptf(:,4) == 9,:);
            inptPTyx = inptf(inptf(:,4) == 10,:);
            inptZ    = [inptZ;inptPTxx;inptPTyy;inptPTxy;inptPTyx];

        end  
        npara = 4;
end
d = inptZ(:,6) + 1i*inptZ(:,7);

Aiqt_Comm = cell(Num_Comm,1);        
for icomm = 1:Num_Comm
    Aiqt_Comm{icomm} = ReadSol_MultRHS(npara*nrec*nfreq_per_comm(icomm),['AIQT_Comm_' num2str(icomm-1)],0,compress);  
end
AIqt = [];        
for icomm = 1:Num_Comm
    AIqt = [AIqt,Aiqt_Comm{icomm}];
end
clear Aiqt_Comm;

JZ = [];
for ifreq=1:nfreq
    
    AiqtEH = cell(npara,1);
    for ipara = 1:npara
        AiqtEH{ipara} = AIqt(:,(ifreq-1)*npara*nrec + (ipara - 1)*nrec + 1:(ifreq-1)*npara*nrec + ipara*nrec);
    end
    
    switch inv_flag
        case 1 % off-diagonal impedance
            
            JEx    = cell(ns,1);
            JEy    = cell(ns,1);
            JHx    = cell(ns,1);
            JHy    = cell(ns,1);
            for isrc = 1:ns
                 Gb = dbdm(omega(ifreq),minv_ind,MeshTop.TE(minv_ind,:),...
                    MeshTop.ed,Ni,EPGauss{ifreq,isrc}(minv_ind,:,:,:),Gauss3DW,JacobDet(minv_ind,:),Edge_ind);
                 Gk = getGk(omega(1),omega(ifreq),Dm3D,MeshTop.TE(minv_ind,:),MeshTop.ed,M(:,(ifreq-1)*ns+isrc),Edge_ind,minv_ind);
                 G1 = (Gb - Gk);

                 JEx{isrc} = (G1.'*AiqtEH{1}).';
                 JEy{isrc} = (G1.'*AiqtEH{2}).';

                 JHx{isrc} = -1./(1i*omega(ifreq)*mu0)*(G1.'*AiqtEH{3}).'; 
                 JHy{isrc} = -1./(1i*omega(ifreq)*mu0)*(G1.'*AiqtEH{4}).';

            end
            clear AiqtEH
            
            intoutf = intoutFEM(intoutFEM(:,6)==freq(ifreq),:);
            indEx1  = find(intoutf(:,5)==1 & intoutf(:,4)==1);
            indEy1  = find(intoutf(:,5)==1 & intoutf(:,4)==2);
            indHx1  = find(intoutf(:,5)==1 & intoutf(:,4)==4);
            indHy1  = find(intoutf(:,5)==1 & intoutf(:,4)==5);

            indEx2  = find(intoutf(:,5)==2 & intoutf(:,4)==1);
            indEy2  = find(intoutf(:,5)==2 & intoutf(:,4)==2);
            indHx2  = find(intoutf(:,5)==2 & intoutf(:,4)==4);
            indHy2  = find(intoutf(:,5)==2 & intoutf(:,4)==5);

            Ex1     = (intoutf(indEx1,7)+intoutf(indEx1,9)) + 1i*(intoutf(indEx1,8)+intoutf(indEx1,10));
            Ey1     = (intoutf(indEy1,7)+intoutf(indEy1,9)) + 1i*(intoutf(indEy1,8)+intoutf(indEy1,10));
            Hx1     = (intoutf(indHx1,7)+intoutf(indHx1,9)) + 1i*(intoutf(indHx1,8)+intoutf(indHx1,10));
            Hy1     = (intoutf(indHy1,7)+intoutf(indHy1,9)) + 1i*(intoutf(indHy1,8)+intoutf(indHy1,10));

            Ex2     = (intoutf(indEx2,7)+intoutf(indEx2,9)) + 1i*(intoutf(indEx2,8)+intoutf(indEx2,10));
            Ey2     = (intoutf(indEy2,7)+intoutf(indEy2,9)) + 1i*(intoutf(indEy2,8)+intoutf(indEy2,10));
            Hx2     = (intoutf(indHx2,7)+intoutf(indHx2,9)) + 1i*(intoutf(indHx2,8)+intoutf(indHx2,10));
            Hy2     = (intoutf(indHy2,7)+intoutf(indHy2,9)) + 1i*(intoutf(indHy2,8)+intoutf(indHy2,10));
            
            JZxy    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hx1.*JEx{2} + Ex2.*JHx{1} - Hx2.*JEx{1} - Ex1.*JHx{2})...
                    -  (Ex2.*Hx1 - Ex1.*Hx2).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2; 

            JZyx    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hy2.*JEy{1} + Ey1.*JHy{2} - Hy1.*JEy{2} - Ey2.*JHy{1})...
                    -  (Ey1.*Hy2 - Ey2.*Hy1).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2;
            
            clear JEx JEy JHx JHy
            
            JZxy    = JZxy.*(sigx(minv_ind)');
            JZyx    = JZyx.*(sigx(minv_ind)'); 
            JZ      = [JZ;JZxy;JZyx];
            
        case 2 % full impedance
            
            JEx    = cell(ns,1);
            JEy    = cell(ns,1);
            JHx    = cell(ns,1);
            JHy    = cell(ns,1);
            for isrc = 1:ns
                 Gb = dbdm(omega(ifreq),minv_ind,MeshTop.TE(minv_ind,:),...
                    MeshTop.ed,Ni,EPGauss{ifreq,isrc}(minv_ind,:,:,:),Gauss3DW,JacobDet(minv_ind,:),Edge_ind);
                 Gk = getGk(omega(1),omega(ifreq),Dm3D,MeshTop.TE(minv_ind,:),MeshTop.ed,M(:,(ifreq-1)*ns+isrc),Edge_ind,minv_ind);
                 G1 = (Gb - Gk);

                 JEx{isrc} = (G1.'*AiqtEH{1}).';
                 JEy{isrc} = (G1.'*AiqtEH{2}).';

                 JHx{isrc} = -1./(1i*omega(ifreq)*mu0)*(G1.'*AiqtEH{3}).'; 
                 JHy{isrc} = -1./(1i*omega(ifreq)*mu0)*(G1.'*AiqtEH{4}).';

            end
            clear AiqtEH
            
            intoutf = intoutFEM(intoutFEM(:,6)==freq(ifreq),:);
            indEx1  = find(intoutf(:,5)==1 & intoutf(:,4)==1);
            indEy1  = find(intoutf(:,5)==1 & intoutf(:,4)==2);
            indHx1  = find(intoutf(:,5)==1 & intoutf(:,4)==4);
            indHy1  = find(intoutf(:,5)==1 & intoutf(:,4)==5);

            indEx2  = find(intoutf(:,5)==2 & intoutf(:,4)==1);
            indEy2  = find(intoutf(:,5)==2 & intoutf(:,4)==2);
            indHx2  = find(intoutf(:,5)==2 & intoutf(:,4)==4);
            indHy2  = find(intoutf(:,5)==2 & intoutf(:,4)==5);

            Ex1     = (intoutf(indEx1,7)+intoutf(indEx1,9)) + 1i*(intoutf(indEx1,8)+intoutf(indEx1,10));
            Ey1     = (intoutf(indEy1,7)+intoutf(indEy1,9)) + 1i*(intoutf(indEy1,8)+intoutf(indEy1,10));
            Hx1     = (intoutf(indHx1,7)+intoutf(indHx1,9)) + 1i*(intoutf(indHx1,8)+intoutf(indHx1,10));
            Hy1     = (intoutf(indHy1,7)+intoutf(indHy1,9)) + 1i*(intoutf(indHy1,8)+intoutf(indHy1,10));

            Ex2     = (intoutf(indEx2,7)+intoutf(indEx2,9)) + 1i*(intoutf(indEx2,8)+intoutf(indEx2,10));
            Ey2     = (intoutf(indEy2,7)+intoutf(indEy2,9)) + 1i*(intoutf(indEy2,8)+intoutf(indEy2,10));
            Hx2     = (intoutf(indHx2,7)+intoutf(indHx2,9)) + 1i*(intoutf(indHx2,8)+intoutf(indHx2,10));
            Hy2     = (intoutf(indHy2,7)+intoutf(indHy2,9)) + 1i*(intoutf(indHy2,8)+intoutf(indHy2,10));
            
            JZxy    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hx1.*JEx{2} + Ex2.*JHx{1} - Hx2.*JEx{1} - Ex1.*JHx{2})...
                    -  (Ex2.*Hx1 - Ex1.*Hx2).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2; 

            JZyx    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hy2.*JEy{1} + Ey1.*JHy{2} - Hy1.*JEy{2} - Ey2.*JHy{1})...
                    -  (Ey1.*Hy2 - Ey2.*Hy1).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2;

            JZxx    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hy2.*JEx{1} + Ex1.*JHy{2} - Hy1.*JEx{2} - Ex2.*JHy{1})...
                    -  (Ex1.*Hy2 - Ex2.*Hy1).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2;

            JZyy    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hx1.*JEy{2} + Ey2.*JHx{1} - Hx2.*JEy{1} - Ey1.*JHx{2})...
                    -  (Ey2.*Hx1 - Ey1.*Hx2).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2;

            clear JEx JEy JHx JHy

            JZxy    = JZxy.*(sigx(minv_ind)');
            JZyx    = JZyx.*(sigx(minv_ind)'); 
            JZxx    = JZxx.*(sigx(minv_ind)');
            JZyy    = JZyy.*(sigx(minv_ind)'); 
            JZ      = [JZ;JZxy;JZyx;JZxx;JZyy];
            
        case 3 % tipper
            
            JHx    = cell(ns,1);
            JHy    = cell(ns,1);
            JHz    = cell(ns,1);
            for isrc = 1:ns
                Gb = dbdm(omega(ifreq),minv_ind,MeshTop.TE(minv_ind,:),...
                    MeshTop.ed,Ni,EPGauss{ifreq,isrc}(minv_ind,:,:,:),Gauss3DW,JacobDet(minv_ind,:),Edge_ind);
                 Gk = getGk(omega(1),omega(ifreq),Dm3D,MeshTop.TE(minv_ind,:),MeshTop.ed,M(:,(ifreq-1)*ns+isrc),Edge_ind,minv_ind);
                 G1 = (Gb - Gk);
                 
                 JHx{isrc} = -1./(1i*omega(ifreq)*mu0)*(G1.'*AiqtEH{1}).'; 
                 JHy{isrc} = -1./(1i*omega(ifreq)*mu0)*(G1.'*AiqtEH{2}).';
                 JHz{isrc} = -1./(1i*omega(ifreq)*mu0)*(G1.'*AiqtEH{3}).';
            end
            clear AiqtEH
            
            intoutf = intoutFEM(intoutFEM(:,6)==freq(ifreq),:);
            indHx1  = find(intoutf(:,5)==1 & intoutf(:,4)==4);
            indHy1  = find(intoutf(:,5)==1 & intoutf(:,4)==5);
            indHz1  = find(intoutf(:,5)==1 & intoutf(:,4)==6);

            indHx2  = find(intoutf(:,5)==2 & intoutf(:,4)==4);
            indHy2  = find(intoutf(:,5)==2 & intoutf(:,4)==5);
            indHz2  = find(intoutf(:,5)==2 & intoutf(:,4)==6);

            Hx1     = (intoutf(indHx1,7)+intoutf(indHx1,9)) + 1i*(intoutf(indHx1,8)+intoutf(indHx1,10));
            Hy1     = (intoutf(indHy1,7)+intoutf(indHy1,9)) + 1i*(intoutf(indHy1,8)+intoutf(indHy1,10));
            Hz1     = (intoutf(indHz1,7)+intoutf(indHz1,9)) + 1i*(intoutf(indHz1,8)+intoutf(indHz1,10));

            Hx2     = (intoutf(indHx2,7)+intoutf(indHx2,9)) + 1i*(intoutf(indHx2,8)+intoutf(indHx2,10));
            Hy2     = (intoutf(indHy2,7)+intoutf(indHy2,9)) + 1i*(intoutf(indHy2,8)+intoutf(indHy2,10));
            Hz2     = (intoutf(indHz2,7)+intoutf(indHz2,9)) + 1i*(intoutf(indHz2,8)+intoutf(indHz2,10));
            
            JTx    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hy2.*JHz{1} + Hz1.*JHy{2} - Hy1.*JHz{2} - Hz2.*JHy{1})...
                    -  (Hz1.*Hy2 - Hz2.*Hy1).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2;
                
            JTy    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hx1.*JHz{2} + Hz2.*JHx{1} - Hx2.*JHz{1} - Hz1.*JHx{2})...
                    -  (Hz2.*Hx1 - Hz1.*Hx2).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2;
            
            clear JHx JHy JHz
            
            JTx    = JTx.*(sigx(minv_ind)');
            JTy    = JTy.*(sigx(minv_ind)');
            JZ     = [JZ;JTx;JTy];
            
        case 4 % off diagonal impedance + tipper
            
            JEx    = cell(ns,1);
            JEy    = cell(ns,1);
            JHx    = cell(ns,1);
            JHy    = cell(ns,1);
            JHz    = cell(ns,1);
            for isrc = 1:ns
                 Gb = dbdm(omega(ifreq),minv_ind,MeshTop.TE(minv_ind,:),...
                    MeshTop.ed,Ni,EPGauss{ifreq,isrc}(minv_ind,:,:,:),Gauss3DW,JacobDet(minv_ind,:),Edge_ind);
                 Gk = getGk(omega(1),omega(ifreq),Dm3D,MeshTop.TE(minv_ind,:),MeshTop.ed,M(:,(ifreq-1)*ns+isrc),Edge_ind,minv_ind);
                 G1 = (Gb - Gk);

                 JEx{isrc} = (G1.'*AiqtEH{1}).';
                 JEy{isrc} = (G1.'*AiqtEH{2}).';

                 JHx{isrc} = -1./(1i*omega(ifreq)*mu0)*(G1.'*AiqtEH{3}).'; 
                 JHy{isrc} = -1./(1i*omega(ifreq)*mu0)*(G1.'*AiqtEH{4}).';
                 JHz{isrc} = -1./(1i*omega(ifreq)*mu0)*(G1.'*AiqtEH{5}).';

            end
            clear AiqtEH
            
            intoutf = intoutFEM(intoutFEM(:,6)==freq(ifreq),:);
            indEx1  = find(intoutf(:,5)==1 & intoutf(:,4)==1);
            indEy1  = find(intoutf(:,5)==1 & intoutf(:,4)==2);
            indHx1  = find(intoutf(:,5)==1 & intoutf(:,4)==4);
            indHy1  = find(intoutf(:,5)==1 & intoutf(:,4)==5);
            indHz1  = find(intoutf(:,5)==1 & intoutf(:,4)==6);

            indEx2  = find(intoutf(:,5)==2 & intoutf(:,4)==1);
            indEy2  = find(intoutf(:,5)==2 & intoutf(:,4)==2);
            indHx2  = find(intoutf(:,5)==2 & intoutf(:,4)==4);
            indHy2  = find(intoutf(:,5)==2 & intoutf(:,4)==5);
            indHz2  = find(intoutf(:,5)==2 & intoutf(:,4)==6);

            Ex1     = (intoutf(indEx1,7)+intoutf(indEx1,9)) + 1i*(intoutf(indEx1,8)+intoutf(indEx1,10));
            Ey1     = (intoutf(indEy1,7)+intoutf(indEy1,9)) + 1i*(intoutf(indEy1,8)+intoutf(indEy1,10));
            Hx1     = (intoutf(indHx1,7)+intoutf(indHx1,9)) + 1i*(intoutf(indHx1,8)+intoutf(indHx1,10));
            Hy1     = (intoutf(indHy1,7)+intoutf(indHy1,9)) + 1i*(intoutf(indHy1,8)+intoutf(indHy1,10));
            Hz1     = (intoutf(indHz1,7)+intoutf(indHz1,9)) + 1i*(intoutf(indHz1,8)+intoutf(indHz1,10));

            Ex2     = (intoutf(indEx2,7)+intoutf(indEx2,9)) + 1i*(intoutf(indEx2,8)+intoutf(indEx2,10));
            Ey2     = (intoutf(indEy2,7)+intoutf(indEy2,9)) + 1i*(intoutf(indEy2,8)+intoutf(indEy2,10));
            Hx2     = (intoutf(indHx2,7)+intoutf(indHx2,9)) + 1i*(intoutf(indHx2,8)+intoutf(indHx2,10));
            Hy2     = (intoutf(indHy2,7)+intoutf(indHy2,9)) + 1i*(intoutf(indHy2,8)+intoutf(indHy2,10));
            Hz2     = (intoutf(indHz2,7)+intoutf(indHz2,9)) + 1i*(intoutf(indHz2,8)+intoutf(indHz2,10));
            
            JZxy    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hx1.*JEx{2} + Ex2.*JHx{1} - Hx2.*JEx{1} - Ex1.*JHx{2})...
                    -  (Ex2.*Hx1 - Ex1.*Hx2).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2; 

            JZyx    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hy2.*JEy{1} + Ey1.*JHy{2} - Hy1.*JEy{2} - Ey2.*JHy{1})...
                    -  (Ey1.*Hy2 - Ey2.*Hy1).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2;
            
            JTx    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hy2.*JHz{1} + Hz1.*JHy{2} - Hy1.*JHz{2} - Hz2.*JHy{1})...
                    -  (Hz1.*Hy2 - Hz2.*Hy1).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2;
                
            JTy    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hx1.*JHz{2} + Hz2.*JHx{1} - Hx2.*JHz{1} - Hz1.*JHx{2})...
                    -  (Hz2.*Hx1 - Hz1.*Hx2).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2;
                
            clear JEx JEy JHx JHy JHz
            
            JZxy    = JZxy.*(sigx(minv_ind)');
            JZyx    = JZyx.*(sigx(minv_ind)'); 
            JTx     = JTx.*(sigx(minv_ind)');
            JTy     = JTy.*(sigx(minv_ind)');
            JZ      = [JZ;JZxy;JZyx;JTx;JTy];
            
        case 5 % full impedance + tipper
            
            JEx    = cell(ns,1);
            JEy    = cell(ns,1);
            JHx    = cell(ns,1);
            JHy    = cell(ns,1);
            JHz    = cell(ns,1);
            for isrc = 1:ns
                 Gb = dbdm(omega(ifreq),minv_ind,MeshTop.TE(minv_ind,:),...
                    MeshTop.ed,Ni,EPGauss{ifreq,isrc}(minv_ind,:,:,:),Gauss3DW,JacobDet(minv_ind,:),Edge_ind);
                 Gk = getGk(omega(1),omega(ifreq),Dm3D,MeshTop.TE(minv_ind,:),MeshTop.ed,M(:,(ifreq-1)*ns+isrc),Edge_ind,minv_ind);
                 G1 = (Gb - Gk);

                 JEx{isrc} = (G1.'*AiqtEH{1}).';
                 JEy{isrc} = (G1.'*AiqtEH{2}).';

                 JHx{isrc} = -1./(1i*omega(ifreq)*mu0)*(G1.'*AiqtEH{3}).'; 
                 JHy{isrc} = -1./(1i*omega(ifreq)*mu0)*(G1.'*AiqtEH{4}).';
                 JHz{isrc} = -1./(1i*omega(ifreq)*mu0)*(G1.'*AiqtEH{5}).';

            end
            clear AiqtEH
            
            intoutf = intoutFEM(intoutFEM(:,6)==freq(ifreq),:);
            indEx1  = find(intoutf(:,5)==1 & intoutf(:,4)==1);
            indEy1  = find(intoutf(:,5)==1 & intoutf(:,4)==2);
            indHx1  = find(intoutf(:,5)==1 & intoutf(:,4)==4);
            indHy1  = find(intoutf(:,5)==1 & intoutf(:,4)==5);
            indHz1  = find(intoutf(:,5)==1 & intoutf(:,4)==6);

            indEx2  = find(intoutf(:,5)==2 & intoutf(:,4)==1);
            indEy2  = find(intoutf(:,5)==2 & intoutf(:,4)==2);
            indHx2  = find(intoutf(:,5)==2 & intoutf(:,4)==4);
            indHy2  = find(intoutf(:,5)==2 & intoutf(:,4)==5);
            indHz2  = find(intoutf(:,5)==2 & intoutf(:,4)==6);

            Ex1     = (intoutf(indEx1,7)+intoutf(indEx1,9)) + 1i*(intoutf(indEx1,8)+intoutf(indEx1,10));
            Ey1     = (intoutf(indEy1,7)+intoutf(indEy1,9)) + 1i*(intoutf(indEy1,8)+intoutf(indEy1,10));
            Hx1     = (intoutf(indHx1,7)+intoutf(indHx1,9)) + 1i*(intoutf(indHx1,8)+intoutf(indHx1,10));
            Hy1     = (intoutf(indHy1,7)+intoutf(indHy1,9)) + 1i*(intoutf(indHy1,8)+intoutf(indHy1,10));
            Hz1     = (intoutf(indHz1,7)+intoutf(indHz1,9)) + 1i*(intoutf(indHz1,8)+intoutf(indHz1,10));

            Ex2     = (intoutf(indEx2,7)+intoutf(indEx2,9)) + 1i*(intoutf(indEx2,8)+intoutf(indEx2,10));
            Ey2     = (intoutf(indEy2,7)+intoutf(indEy2,9)) + 1i*(intoutf(indEy2,8)+intoutf(indEy2,10));
            Hx2     = (intoutf(indHx2,7)+intoutf(indHx2,9)) + 1i*(intoutf(indHx2,8)+intoutf(indHx2,10));
            Hy2     = (intoutf(indHy2,7)+intoutf(indHy2,9)) + 1i*(intoutf(indHy2,8)+intoutf(indHy2,10));
            Hz2     = (intoutf(indHz2,7)+intoutf(indHz2,9)) + 1i*(intoutf(indHz2,8)+intoutf(indHz2,10));
            
            JZxy    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hx1.*JEx{2} + Ex2.*JHx{1} - Hx2.*JEx{1} - Ex1.*JHx{2})...
                    -  (Ex2.*Hx1 - Ex1.*Hx2).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2; 

            JZyx    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hy2.*JEy{1} + Ey1.*JHy{2} - Hy1.*JEy{2} - Ey2.*JHy{1})...
                    -  (Ey1.*Hy2 - Ey2.*Hy1).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2;

            JZxx    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hy2.*JEx{1} + Ex1.*JHy{2} - Hy1.*JEx{2} - Ex2.*JHy{1})...
                    -  (Ex1.*Hy2 - Ex2.*Hy1).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2;

            JZyy    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hx1.*JEy{2} + Ey2.*JHx{1} - Hx2.*JEy{1} - Ey1.*JHx{2})...
                    -  (Ey2.*Hx1 - Ey1.*Hx2).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2;
            
            JTx    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hy2.*JHz{1} + Hz1.*JHy{2} - Hy1.*JHz{2} - Hz2.*JHy{1})...
                    -  (Hz1.*Hy2 - Hz2.*Hy1).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2;
                
            JTy    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hx1.*JHz{2} + Hz2.*JHx{1} - Hx2.*JHz{1} - Hz1.*JHx{2})...
                    -  (Hz2.*Hx1 - Hz1.*Hx2).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2;
                
            clear JEx JEy JHx JHy JHz
            
            JZxy    = JZxy.*(sigx(minv_ind)');
            JZyx    = JZyx.*(sigx(minv_ind)');
            JZxx    = JZxx.*(sigx(minv_ind)');
            JZyy    = JZyy.*(sigx(minv_ind)');
            JTx     = JTx.*(sigx(minv_ind)');
            JTy     = JTy.*(sigx(minv_ind)');
            JZ      = [JZ;JZxy;JZyx;JZxx;JZyy;JTx;JTy];
            
        case 6 % phase tensor
            
            JEx    = cell(ns,1);
            JEy    = cell(ns,1);
            JHx    = cell(ns,1);
            JHy    = cell(ns,1);
            for isrc = 1:ns
                 Gb = dbdm(omega(ifreq),minv_ind,MeshTop.TE(minv_ind,:),...
                    MeshTop.ed,Ni,EPGauss{ifreq,isrc}(minv_ind,:,:,:),Gauss3DW,JacobDet(minv_ind,:),Edge_ind);
                 Gk = getGk(omega(1),omega(ifreq),Dm3D,MeshTop.TE(minv_ind,:),MeshTop.ed,M(:,(ifreq-1)*ns+isrc),Edge_ind,minv_ind);
                 G1 = (Gb - Gk);

                 JEx{isrc} = (G1.'*AiqtEH{1}).';
                 JEy{isrc} = (G1.'*AiqtEH{2}).';

                 JHx{isrc} = -1./(1i*omega(ifreq)*mu0)*(G1.'*AiqtEH{3}).'; 
                 JHy{isrc} = -1./(1i*omega(ifreq)*mu0)*(G1.'*AiqtEH{4}).';

            end
            clear AiqtEH
            
            intoutf = intoutFEM(intoutFEM(:,6)==freq(ifreq),:);
            indEx1  = find(intoutf(:,5)==1 & intoutf(:,4)==1);
            indEy1  = find(intoutf(:,5)==1 & intoutf(:,4)==2);
            indHx1  = find(intoutf(:,5)==1 & intoutf(:,4)==4);
            indHy1  = find(intoutf(:,5)==1 & intoutf(:,4)==5);

            indEx2  = find(intoutf(:,5)==2 & intoutf(:,4)==1);
            indEy2  = find(intoutf(:,5)==2 & intoutf(:,4)==2);
            indHx2  = find(intoutf(:,5)==2 & intoutf(:,4)==4);
            indHy2  = find(intoutf(:,5)==2 & intoutf(:,4)==5);

            Ex1     = (intoutf(indEx1,7)+intoutf(indEx1,9)) + 1i*(intoutf(indEx1,8)+intoutf(indEx1,10));
            Ey1     = (intoutf(indEy1,7)+intoutf(indEy1,9)) + 1i*(intoutf(indEy1,8)+intoutf(indEy1,10));
            Hx1     = (intoutf(indHx1,7)+intoutf(indHx1,9)) + 1i*(intoutf(indHx1,8)+intoutf(indHx1,10));
            Hy1     = (intoutf(indHy1,7)+intoutf(indHy1,9)) + 1i*(intoutf(indHy1,8)+intoutf(indHy1,10));

            Ex2     = (intoutf(indEx2,7)+intoutf(indEx2,9)) + 1i*(intoutf(indEx2,8)+intoutf(indEx2,10));
            Ey2     = (intoutf(indEy2,7)+intoutf(indEy2,9)) + 1i*(intoutf(indEy2,8)+intoutf(indEy2,10));
            Hx2     = (intoutf(indHx2,7)+intoutf(indHx2,9)) + 1i*(intoutf(indHx2,8)+intoutf(indHx2,10));
            Hy2     = (intoutf(indHy2,7)+intoutf(indHy2,9)) + 1i*(intoutf(indHy2,8)+intoutf(indHy2,10));
    
            Zxx      = (Ex1.*Hy2 - Ex2.*Hy1)./(Hx1.*Hy2-Hx2.*Hy1);
            Zxy      = (Ex2.*Hx1 - Ex1.*Hx2)./(Hx1.*Hy2 - Hx2.*Hy1);
            Zyx      = (Ey1.*Hy2 - Ey2.*Hy1)./(Hx1.*Hy2-Hx2.*Hy1);
            Zyy      = (Ey2.*Hx1 - Ey1.*Hx2)./(Hx1.*Hy2-Hx2.*Hy1);

            JZxy    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hx1.*JEx{2} + Ex2.*JHx{1} - Hx2.*JEx{1} - Ex1.*JHx{2})...
                    -  (Ex2.*Hx1 - Ex1.*Hx2).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2; 

            JZyx    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hy2.*JEy{1} + Ey1.*JHy{2} - Hy1.*JEy{2} - Ey2.*JHy{1})...
                    -  (Ey1.*Hy2 - Ey2.*Hy1).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2;

            JZxx    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hy2.*JEx{1} + Ex1.*JHy{2} - Hy1.*JEx{2} - Ex2.*JHy{1})...
                    -  (Ex1.*Hy2 - Ex2.*Hy1).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2;

            JZyy    = ((Hx1.*Hy2 - Hx2.*Hy1).*(Hx1.*JEy{2} + Ey2.*JHx{1} - Hx2.*JEy{1} - Ey1.*JHx{2})...
                    -  (Ey2.*Hx1 - Ey1.*Hx2).*(Hy2.*JHx{1} + Hx1.*JHy{2} - Hy1.*JHx{2} - Hx2.*JHy{1}))./(Hx1.*Hy2 - Hx2.*Hy1).^2;
            
            clear JEx JEy JHx JHy
            
            JZxy    = JZxy.*(sigx(minv_ind)');
            JZyx    = JZyx.*(sigx(minv_ind)'); 
            JZxx    = JZxx.*(sigx(minv_ind)');
            JZyy    = JZyy.*(sigx(minv_ind)');

            JPTxx   = ((real(Zxx).*real(Zyy) - real(Zxy).*real(Zyx)).*(real(JZyy).*imag(Zxx) + real(Zyy).*imag(JZxx) ...
                        - real(JZxy).*imag(Zyx) - real(Zxy).*imag(JZyx)) - (real(Zyy).*imag(Zxx) - real(Zxy).*imag(Zyx)).* ...
                        (real(JZxx).*real(Zyy) + real(Zxx).*real(JZyy) - real(JZxy).*real(Zyx) - real(Zxy).*real(JZyx)))./(real(Zxx).*real(Zyy) - real(Zxy).*real(Zyx)).^2;
            
            JPTyy   = ((real(Zxx).*real(Zyy) - real(Zxy).*real(Zyx)).*(-real(JZyx).*imag(Zxy) - real(Zyx).*imag(JZxy) ...
                        + real(JZxx).*imag(Zyy) + real(Zxx).*imag(JZyy)) - (-real(Zyx).*imag(Zxy) + real(Zxx).*imag(Zyy)).* ...
                        (real(JZxx).*real(Zyy) + real(Zxx).*real(JZyy) - real(JZxy).*real(Zyx) - real(Zxy).*real(JZyx)))./(real(Zxx).*real(Zyy) - real(Zxy).*real(Zyx)).^2;

            JPTxy   = ((real(Zxx).*real(Zyy) - real(Zxy).*real(Zyx)).*(real(JZyy).*imag(Zxy) + real(Zyy).*imag(JZxy) ...
                    - real(JZxy).*imag(Zyy) - real(Zxy).*imag(JZyy)) - (real(Zyy).*imag(Zxy) - real(Zxy).*imag(Zyy)).* ...
                    (real(JZxx).*real(Zyy) + real(Zxx).*real(JZyy) - real(JZxy).*real(Zyx) - real(Zxy).*real(JZyx)))./(real(Zxx).*real(Zyy) - real(Zxy).*real(Zyx)).^2;

            JPTyx   = ((real(Zxx).*real(Zyy) - real(Zxy).*real(Zyx)).*(-real(JZyx).*imag(Zxx) - real(Zyx).*imag(JZxx) ...
                        + real(JZxx).*imag(Zyx) + real(Zxx).*imag(JZyx)) - (-real(Zyx).*imag(Zxx) + real(Zxx).*imag(Zyx)).* ...
                        (real(JZxx).*real(Zyy) + real(Zxx).*real(JZyy) - real(JZxy).*real(Zyx) - real(Zxy).*real(JZyx)))./(real(Zxx).*real(Zyy) - real(Zxy).*real(Zyx)).^2;
            
            clear JZxy JZyx JZxx JZyy
            
            JZ      = [JZ;JPTxx;JPTyy;JPTxy;JPTyx];
            
    end
    
end
