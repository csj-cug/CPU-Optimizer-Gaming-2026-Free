% function CSEM3DFOR_Hex_MumpsPardiso
%%update: direct solver is added
%%The curl of Ni is saved to be used for applying Farady's law to get
%%magnetic field. No numerical differential is required to get magnetic
%%field from electric field in such case. The electric field is defined at
%%the center of each cell so far. (see function getTE_Hxyz_CellCenter),
%%magnetic field at the receivers will be interpolated from magnetic field
%%in the center of each cell

%  addpath /public/home/caihongzhu/zhuhao/MT3D/MT3D_structured_Inv/green3d/double_aug13_intel_omp/bin
%  addpath /public/home/caihongzhu/zhuhao/MT3D/MT3D_structured_Inv/hdf5/matlab
 addpath /public/home/caihongzhu/zhuhao/Joint/Joint_structured_Inv/green3d/double_aug13_intel_omp/bin
  addpath /public/home/caihongzhu/zhuhao/Joint/Joint_structured_Inv/hdf5/matlab
%addpath /public/home/caihongzhu/zhuhao/Joint/Joint_structured_Inv/input2
tic;
global  mu0;
global  eps0;
global  sig_air;
Input;

save anomalousBnd xa ya za

if ~iscell(par)       
    error('par should be a cell structure');
end
ns = length(par);              
load recpar.dat;

if length(unique(recpar(:,5))) ~= ns
    error('Numer of source in recpar.dat does not conform with the Input file');
end

freq    = unique(recpar(:,6));
nfreq   = length(freq);
omega   = 2*pi*freq;
mu0     = 4*pi*1e-7;
eps0    = 8.854*1e-12;
sig_air = 1e-9;

fprintf('%s \n','Stage 1: Prepare for nodes and connectivity matrix');
[X, Y, Z] = ndgrid(x,y,z);
if exist('top.mat','file')
    fprintf('Load topography file \n');
    load top.mat
else
    fprintf('top file is not found, give up this option \n');
    [Xtop,Ytop,Ztop] = XYZTop(x,y);     %Load seafloor bathymetry
end
MeshReg = connectivity_brick_ByNode(X,Y,Z);

if length(unique(Ztop(:))) == 1 && unique(Ztop(:))==ShiftZ
    warning('Input topography is flat');   
    MeshTop = MeshReg;
    MeshTop.Ztop_interp = ShiftZ*ones(size(X(:,:,1)));
    MeshTop.Z_new = Z;
else
    MeshTop = connectivity_brick_Topo_New(X,Y,Z,Xtop,Ytop,Ztop,BoundUp,BoundLow,ShiftZ,AspectChangeRange);
end

bnd                                                       = getedge_bnd(MeshTop.ed_center); %logical
[Gauss, GaussW, Gauss3Dx, Gauss3Dy, Gauss3Dz, Gauss3DW]   = GetGaussian(GaussOrder);
L_cell                                                    = get_NodeBasis_OneCell(Gauss3Dx,Gauss3Dy,Gauss3Dz);     
[Jacob, JacobInv, JacobDet]                               = GetJacobian(MeshTop.p,MeshTop.t,Gauss);

[sig0x, sig0y, sig0z, sigx, sigy, sigz] = SetSig_New(MeshTop.t_xyz_center,MeshReg.t_xyz_center,xa,ya,za,sigax,sigay,sigaz,...
    hl,sl,al,loadSig,AnoGridReg,viewflag);
save conductivity_true sig0x sig0y sig0z sigx sigy sigz

epxnode = cell(nfreq,ns);epynode = cell(nfreq,ns);epznode = cell(nfreq,ns);

for ifreq = 1:nfreq
    for isrc = 1:ns
        [epxnode{ifreq,isrc},epynode{ifreq,isrc},epznode{ifreq,isrc}] = getprimaryGreen_Par(freq(ifreq),par{isrc},MeshTop.p,sl,hl,al,viewflag);
    end
end

EPGauss = cell(nfreq,ns);

 for ifreq = 1:nfreq
    for isrc = 1:ns
        EPGauss{ifreq,isrc} = GetXYZGauss(MeshTop.p,MeshTop.t,L_cell,epxnode{ifreq,isrc},epynode{ifreq,isrc},epznode{ifreq,isrc});
    end
end

Ve = getvolume(JacobDet,Gauss3DW);

fprintf('%s \n','Stage 2: Produce local/global matrix and vector');
fprintf('%s \n','Only the matrix for the first frequency will be produced');

[K,~,Ni,CurlNi] = getK_Brick_Edge_Jacobian(omega(1),sigx,sigy,sigz,MeshTop.TE,MeshTop.TE_length,MeshTop.ed,Jacob,JacobInv,JacobDet,Gauss,Gauss3DW);
[rec,ind_rec_in_recpar,ind_recpar_in_rec] = unique(recpar(:,1:3),'rows','stable');
nrec=size(rec,1);
[Q_Ex, Q_Ey, Q_Ez,Q_Hx,Q_Hy,Q_Hz] = GetInterpMatrix_Rec(rec,Ve,MeshReg.t_xyz_mm,Ni,CurlNi,MeshReg.ed_length,MeshReg.TE); %Field Interpolation matrix

b = cell(nfreq,ns);

for ifreq = 1:nfreq
    for isrc = 1:ns
        b{ifreq,isrc} = getB_edge_New(omega(ifreq),sig0x,sig0y,sig0z,sigx,sigy,sigz,MeshTop.TE,MeshTop.ed,Ni,EPGauss{ifreq,isrc},Gauss3DW,JacobDet);
    end
end

fprintf('%s \n','Stage 3: Add boundary condition');
Bdirc         = getbnd(bnd);                                              % Get boundary condition matrix
[Kglob,bglob] = DiricBndNew_MultSrc_MultFreq(K,b,Bdirc,freq,nfreq,ns);    % Apply Dirichlet boundary condition

fprintf('%s \n','Stage 4: Solve the linear equations');
[freq_per_comm, nfreq_per_comm] = PardisoInputGeneralMultRHSMultFreq_Split(Kglob,bglob,1,symFlag,fileformat,compress,freq,Num_Comm);

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

res = zeros(ns*nfreq,1);
m   = [];

for ifreq = 1:nfreq
    for isrc = 1:ns
        res((ifreq-1)*ns+isrc) = norm((real(Kglob)+...
            (freq(ifreq)/freq(1))*1i*imag(Kglob))*M(:,(ifreq-1)*ns+isrc)-bglob{ifreq,isrc})/norm(bglob{ifreq,isrc});
        %       m{isrc}   = M(:,isrc);
        fprintf('%s %d\n',['Normalized residual from ' SOLVERSTR 'for frequency' num2str(freq(ifreq)) ' Hz source ' num2str(isrc) ' is'], res((ifreq-1)*ns+isrc));
    end
end

%Get electromagnetic fields
intoutFEM = [];

for ifreq = 1:nfreq
    for isrc = 1:ns
        intoutFEMtmp = getfield_new(M(:,(ifreq-1)*ns+isrc),Q_Ex,Q_Ey,Q_Ez,Q_Hx,Q_Hy,Q_Hz,recpar,hl,sl,al,par{isrc},freq(ifreq),omega(ifreq),isrc,Xtop,Ytop,Ztop,MeshTop.p);
        intoutFEM = [intoutFEM;intoutFEMtmp];
    end
end

%Get full impedances
[inptdata,inptdataB] = GetImpedanceFull(intoutFEM);

%Add noise
if exist('noise_flag','var')
    if noise_flag == 1
        fprintf('Add %g %g Gauss noise to impedance and tipper data, respectively.\n',noise_lel_Z,noise_lel_T);
        inptdata = AddGaussNoise(inptdata,noise_lel_Z,noise_lel_T);
    elseif noise_flag == 0
        fprintf('no noise is added to synthetic data.\n');
    end
end

%Compute and add standard deviation
[inptdata, inptdataB] = getstandarddeviate(inptdata, inptdataB,errZ_off,errZ_on);

save inptdata.dat inptdata -ascii;
save inptdataB.dat inptdataB -ascii;
