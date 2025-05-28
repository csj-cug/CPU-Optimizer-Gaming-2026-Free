function d = FEM_forward(freq,ns,recpar,par,omega,np,npo,MeshTop,bnd,hl,sl,al,sig0x,sig0y,sig0z,sigx,sigy,...
           sigz,Jacob,JacobInv,JacobDet,Gauss3DW,Gauss,EPGauss,compress,symFlag,fileformat,Num_Comm,Q_Ex,Q_Ey,Q_Ez,Q_Hx,Q_Hy,Q_Hz,mpirun,Xtop,Ytop,Ztop,inv_flag)
       
nfreq=length(freq);
[K, ~, Ni, ~] = getK_Brick_Edge_Jacobian(omega(1),sigx,sigy,sigz,MeshTop.TE,MeshTop.TE_length,MeshTop.ed,Jacob,JacobInv,JacobDet,Gauss,Gauss3DW);
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

mpibin = '/public/home/caihongzhu/zhuhao/Joint/Joint_structured_Inv/pardiso_scratch_MpiSplit_small/intel_intelmpi_lp64_intel64_so';

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
    case 3
        inptZ = [];
        for ifreq = 1:nfreq
            inptf   = inptdata(inptdata(:,5) == freq(ifreq),:);
            inptTx  = inptf(inptf(:,4) == 5,:);
            inptTy  = inptf(inptf(:,4) == 6,:);
            inptZ   = [inptZ;inptTx;inptTy];
        end   
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
end
d        = inptZ(:,6) + 1i*inptZ(:,7);
