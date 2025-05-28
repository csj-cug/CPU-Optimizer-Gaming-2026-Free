%% 3-D magnetotelluric isotropic inversion 
%% Explicit sensitivity matrix
%% Gauss - Newton method
%% smothness-maximum inversion

addpath /public/home/caihongzhu/zhuhao/Joint/Joint_structured_Inv/green3d/double_aug13_intel_omp/bin
addpath /public/home/caihongzhu/zhuhao/Joint/Joint_structured_Inv/hdf5/matlab
%addpath /public/home/caihongzhu/zhuhao/Joint/Joint_structured_Inv/model1/joint_wuxing

% 
clear;
tic;

fprintf('%s \n','Smothness-maximum isotropic inversion based on explicit sensitivity matrix is adopted!');

global  mu0;
global  eps0;
global  sig_air;
InputInv;

if ~iscell(par)
    error('par should be a cell structure');
end
ns = length(par);

load recpar.dat;
load otimes.mat;

[vgrids_ref,vgrids_true] = read_VgridApr;
vgrids  = vgrids_ref;
Creat_New_Vgrids(vgrids);

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

x1 = x;
y1 = y;
z1 = z;

if exist('top.mat','file')
    fprintf('Load topography file \n');
    load top.mat
else
    fprintf('top file is not found, give up this option \n');
    [Xtop,Ytop,Ztop]  = XYZTop(x,y);     %Load seafloor bathymetry
end

MeshReg = connectivity_brick_ByNode(X,Y,Z);

t_xyz_center = MeshReg.t_xyz_center;
minv_ind = find(t_xyz_center(:,1)>=invBnd(1) & t_xyz_center(:,1)<=invBnd(2) & t_xyz_center(:,2)>=invBnd(3)...
                 & t_xyz_center(:,2)<=invBnd(4) & t_xyz_center(:,3)>=invBnd(5) & t_xyz_center(:,3)<=invBnd(6));

if length(unique(Ztop(:))) == 1 && unique(Ztop(:))==ShiftZ
    warning('Input topography is flat');
    MeshTop = MeshReg;
    MeshTop.Ztop_interp = ShiftZ*ones(size(X(:,:,1)));
    MeshTop.Z_new = Z;
else
    MeshTop = connectivity_brick_Topo_New(X,Y,Z,Xtop,Ytop,Ztop,BoundUp,BoundLow,ShiftZ,AspectChangeRange);
end

bnd                                                     = getedge_bnd(MeshTop.ed_center);
[Gauss, GaussW, Gauss3Dx, Gauss3Dy, Gauss3Dz, Gauss3DW] = GetGaussian(GaussOrder);
L_cell                                                  = get_NodeBasis_OneCell(Gauss3Dx,Gauss3Dy,Gauss3Dz);      %get node basis function for one cell at the Gauss point
[Jacob, JacobInv, JacobDet]                             = GetJacobian(MeshTop.p,MeshTop.t,Gauss);
[sig0x, sig0y, sig0z, sigx, sigy, sigz]                 = SetSig_New(MeshTop.t_xyz_center,MeshReg.t_xyz_center,...
    xa,ya,za,sigax,sigay,sigaz,hl,sl,al,loadSig,AnoGridReg,viewflag);

epxnode = cell(nfreq,ns);epynode = cell(nfreq,ns);epznode = cell(nfreq,ns);%modofication: define primary field on node instead of edge center

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
Ve                                        = getvolume(JacobDet,Gauss3DW);
[rec,ind_rec_in_recpar,ind_recpar_in_rec] = unique(recpar(:,1:3),'rows','stable');
nrec=size(rec,1);
[~, ~, Ni, CurlNi] = getK_Brick_Edge_Jacobian(omega(1),sigx,sigy,sigz,...
    MeshTop.TE,MeshTop.TE_length,MeshTop.ed,Jacob,JacobInv,JacobDet,Gauss,Gauss3DW);
%Field Interpolation matrix
[Q_Ex, Q_Ey, Q_Ez,Q_Hx,Q_Hy,Q_Hz,Nre] = GetInterpMatrix_RecInv(rec,Ve,MeshReg.t_xyz_mm,Ni,CurlNi,MeshReg.ed_length,MeshReg.TE,inv_flag); 

fid=fopen('m','w');
fprintf(fid,'%d \n',Nre);
fclose(fid);

Edge_ind = unique(MeshTop.TE(minv_ind,:));             
NEd      = length(Edge_ind);             
fid=fopen('ned','w');
fprintf(fid,'%d \n',NEd);
fclose(fid);
fid=fopen('Edge_ind','w');
fprintf(fid,'%d \n',Edge_ind);
fclose(fid);

Nm=length(minv_ind);
if isempty(Nm)
    error('no cell is selected');
end
ned = size(MeshTop.ed,1);%number of edges
Ne  = size(MeshTop.TE,1);

%load observation data
load inptdata.dat

switch inv_flag
    case 1
        fprintf('Off-diagonal impedance data is used in inversion \n');
        inptZ = [];
        for ifreq = 1:nfreq
            inptf   = inptdata(inptdata(:,5) == freq(ifreq),:);
            inptZxy = inptf(inptf(:,4) == 1,:);
            inptZyx = inptf(inptf(:,4) == 2,:);
            inptZ   = [inptZ;inptZxy;inptZyx];
        end
    case 2
        fprintf('Full impedance data is used in inversion \n');
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
        fprintf('Tipper data is used in inversion \n');
        inptZ = [];
        for ifreq = 1:nfreq
            inptf   = inptdata(inptdata(:,5) == freq(ifreq),:);
            inptTx  = inptf(inptf(:,4) == 5,:);
            inptTy  = inptf(inptf(:,4) == 6,:);
            inptZ   = [inptZ;inptTx;inptTy];
        end   
    case 4
        fprintf('Both off-diagonal impedance and tipper data is used in inversion \n');
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
        fprintf('Both full diagonal impedance and tipper data is used in inversion \n');
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
        fprintf('Full phase tensor data is used in inversion \n');
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

dobs     = inptZ(:,6) + 1i*inptZ(:,7);
dobs1    = otimes;
wd       = 1./inptZ(:,8);
recdata  = inptZ(:,1:5);
data_all = zeros(size(dobs,1),7,CN+1);
Nd       = size(dobs1,1);
Wd       = spdiags(wd,0,length(wd),length(wd));
Wd1      = GetWd1(dobs1,Nd);
%normalization
  idx      = max(Wd*dobs)*5;
  idx1     = max(Wd1*dobs1)+20;
  dobs     = real(dobs)/real(idx)+1i*imag(dobs)/imag(idx);
  dobs1    = dobs1/idx1;

if (Num_Comm > nfreq)
    error('Num of MPI communicators cannot be larger than number of frequency');
end

p            = MeshReg.p;
t            = MeshReg.t;
t_xyz_center = MeshReg.t_xyz_center;
vol          = Ve(minv_ind);
vol_norm     = sqrt(vol./max(vol));
R            = spdiags(vol_norm,0,Nm,Nm);
t_xyz_center_inv = t_xyz_center(minv_ind,:);
if exist('Rough_N_Neighbour','var')
    fprintf(['Use ' num2str(Rough_N_Neighbour) ' neighbours to calculate roughness matrix \n']);
    L = GetRoughenMatNew(p,t,minv_ind,t_xyz_center,vol,Rough_N_Neighbour);
else  
    L = GetRoughenMatNew(p,t,minv_ind,t_xyz_center,vol,6);   
end
L = R*L; %incorporate volume effect in to ocam

%Inversion parameter
misfit  = zeros(1,CN+1);
misfit1 = zeros(1,CN+1);
apha    = zeros(1,CN+1);
apha1   = zeros(1,CN+1);
beta    = zeros(1,CN+1);
betax   = zeros(1,CN+1);
betay   = zeros(1,CN+1);
betaz   = zeros(1,CN+1);
dm      = zeros(2*Nm,1);
dm1     = zeros(Nm,1);
dm2     = zeros(Nm,1);
dminv   = zeros(Nm,CN+1);
dminv1  = zeros(Nm,CN+1);
Minv    = zeros(Nm,CN+1);
Minv1   = zeros(Nm,CN+1);

m      = m2Logm(log(sigx(minv_ind)),mmin,mmax)/pp;
m1      = vgrids;
mpre    = m;
mapr    = m;
m1apr   = m1;
m1Log   = m2Logm(m1,m1min,m1max);  
m1aprLog = m2Logm(m1apr,m1min,m1max);
I       = speye(Nm,Nm);
SG      = zeros(1,CN+1);
SGx     = zeros(1,CN+1);
SGy     = zeros(1,CN+1);
SGz     = zeros(1,CN+1);
SWm     = zeros(1,CN+1);
SWm1    = zeros(1,CN+1);
Sd      = zeros(1,CN+1);
Sd1     = zeros(1,CN+1);

if exist('m_initial_flag','var')
    if m_initial_flag==0
        fprintf('User specified that there is no initial model \n');
    else
        fprintf('User specified to use initial model \n');
        if exist('m_initial.mat','file')
            fprintf('Initial model file found \n');
            load m_initial;
            m    = m_initial;
            mpre = m_initial;
            mapr = m_initial;
        else
            fprintf('Initial model file is not found, give up this option \n');
        end
    end
else
    fprintf('m_initial_flag is not found; do not consider initial model \n');
end


     
JacobUpdate = 1;
[Lx,Ly,Lz] = GetLaplance(Nmx,Nmy,Nmz,dmx,dmy,dmz);


for n=1:1+CN
    
    
      sigx(minv_ind) = exp(Logm2m(m,mmin,mmax,pp));
      sigy(minv_ind) = exp(Logm2m(m,mmin,mmax,pp));
      sigz(minv_ind) = exp(Logm2m(m,mmin,mmax,pp));

    if JacobUpdate
        [dpr,J]= FEM_forwardJ(nrec,freq,ns,omega,recpar,par,np,npo,MeshTop,bnd,hl,sl,al,sig0x,sig0y,sig0z,sigx,sigy,...
            sigz,Jacob,JacobInv,JacobDet,Gauss3DW,Gauss,EPGauss,compress,symFlag,fileformat,Num_Comm,Q_Ex,Q_Ey,Q_Ez,Q_Hx,Q_Hy,Q_Hz,mpirun,Edge_ind,minv_ind,Xtop,Ytop,Ztop,inv_flag); 
    else
        dpr  = FEM_forward(freq,ns,recpar,par,omega,np,npo,MeshTop,bnd,hl,sl,al,sig0x,sig0y,sig0z,sigx,sigy,...
            sigz,Jacob,JacobInv,JacobDet,Gauss3DW,Gauss,EPGauss,compress,symFlag,fileformat,Num_Comm,Q_Ex,Q_Ey,Q_Ez,Q_Hx,Q_Hy,Q_Hz,mpirun,Xtop,Ytop,Ztop,inv_flag);
    end
    

    %Verificate sensistivity matrix
%      VerfJnew(J,dpr,freq,ns,recpar,par,omega,np,npo,MeshTop,bnd,hl,sl,al,sig0x,sig0y,sig0z,sigx,sigy,sigz,...
%            Jacob,JacobInv,JacobDet,Gauss3DW,Gauss,EPGauss,compress,symFlag,fileformat,Num_Comm,Q_Ex,Q_Ey,Q_Ez,Q_Hx,Q_Hy,Q_Hz,mpirun,xa,ya,za,minv_ind,Xtop,Ytop,Ztop,inv_flag);
    
    fm3d;
    A_seis     = frechet(Nd,Nm,Nmx,Nmy,Nmz);
    A_seis     = sparse(A_seis);
    mtimes     = read_mtimes;

    J         = real(J)/real(idx)+1i*imag(J)/imag(idx);
    dpr       = real(dpr)/real(idx)+1i*imag(dpr)/imag(idx);
    A_seis    = A_seis/idx1;
    mtimes    = mtimes/idx1;
  
 
    data       = [recdata real(dpr) imag(dpr)];
    dpr1       = dpr;

    rn         = dpr1-dobs; 
    rn1        = mtimes-dobs1;
    misfit(n)  = norm(Wd*rn)/norm(Wd*dobs)
    misfit1(n) = norm(Wd1*rn1)/norm(Wd1*dobs1)
    m1         = flip(reshape(m1,Nmx,Nmy,Nmz),3);  
    m1         = m1(:);                           
    
    dminv(:,n)      = dm1;
    dminv1(:,n)     = dm2;
    Minv(:,n)       = m;
    Minv1(:,n)      = m1;
    data_all(:,:,n) = data;
    
    Lxm1            = Lx*m;
    Lxm2            = Lx*m1;
    Lym1            = Ly*m;
    Lym2            = Ly*m1;
    Lzm1            = Lz*m;
    Lzm2            = Lz*m1;
    Lxyzm1          = (Lx+Ly+Lz)*m;
    Lxyzm2          = (Lx+Ly+Lz)*m1;
    SWm(n)          = norm(L*m);
    SWm1(n)         = norm(L*m1);
    Sd(n)           = norm(Wd*rn);
    Sd1(n)          = norm(Wd1*rn1);
    
    % save inversion result
    save InvResult m m1 dm dm1 dminv dminv1 Minv Minv1 misfit misfit1 apha apha1 x1 y1 z1 minv_ind invBnd t_xyz_center Lxm1 Lxm2 Lym1  Lym2 Lzm1 Lzm2 Lxyzm1 Lxyzm2 SWm SWm1 Sd Sd1 data_all
    
    if (misfit(n)<=final_misfit && misfit1(n)<=final_misfit) || n==CN+1
        fprintf('The convergence condition is reached, exit!\n');
        return;
    end
    if exist('mapr_update','var')
       if  mapr_update==2 && n>=2
           if(norm(misfit(n)-misfit(n-1))/norm(misfit(n-1)) >=0.05)
               mapr= m;
           end
       end
    end
    

    m11 = Logm2m(m,mmin,mmax,pp);
    p1  = spdiags(pp*(mmax-m11).*(m11-mmin)  / (mmax-mmin),0,Nm,Nm);
    p2  = spdiags((m1max-m1).*(m1-m1min) / (m1max-m1min),0,Nm,Nm);
                     
     J         = J*p1;                
     A_seis    = A_seis*p2;          
              
    x   = rand(Nm,1);
    x_1 = rand(Nm,1);
    wj  = Wd*J*x;
    wj1 = Wd1*A_seis*x_1;
    Ax  = J'*Wd*wj;
    Ax1 = A_seis'*Wd1*wj1;
    ao  = norm(Ax)/norm(L'*(L*x));
    ao1 = norm(Ax1)/norm(L'*(L*x_1));
    
    if qq_scale_flag ==1
        apha(n+1) = qq*ao;
        apha1(n+1) = qq*ao1;
    elseif qq_scale_flag == 2
        apha(n+1) = qq*ao/n;
        apha1(n+1) = qq*ao1/n;
    elseif qq_scale_flag == 3
        apha(n+1) = qq*ao/(n^2);
        apha1(n+1) = qq*ao1/(n^2);
    else
        error('No such option for qq_scale_flag \n');
    end
    
    SG(n+1)                       = GetSGnew(m,m1Log,I);
    SGx(n+1)                      = GetSGnew(m,m1Log,Lx);
    SGy(n+1)                      = GetSGnew(m,m1Log,Ly);
    SGz(n+1)                      = GetSGnew(m,m1Log,Lz);
    SGD(n+1)                      = SGx(n+1)+SGy(n+1)+SGz(n+1);
    
    
    switch GramianType
        case 1
             DG           = GetDGnew(m,m1Log,I);
            % Calculate regularization parameter
             beta(n+1)    = GetJointPar_Ap(J,A_seis,Wd,Wd1,DG,qbeta,cbeta,n);
             WG_left      = beta(n+1)*DG;
             WG_right     = beta(n+1)*SG(n+1);
        case 2
             DGx                         = GetDGnew(m,m1Log,Lx);
             DGy                         = GetDGnew(m,m1Log,Ly);
             DGz                         = GetDGnew(m,m1Log,Lz);
%             % Calculate regularization parameter
             betax(n+1)                  = GetJointPar_Ap(J,A_seis,Wd,Wd1,DGx,qbeta,cbeta,n);
             betay(n+1)                  = GetJointPar_Ap(J,A_seis,Wd,Wd1,DGy,qbeta,cbeta,n);
             betaz(n+1)                  = GetJointPar_Ap(J,A_seis,Wd,Wd1,DGz,qbeta,cbeta,n);
             WG_left                     = [ betax(n+1)*DGx; betay(n+1)*DGy; betaz(n+1)*DGz ];
             WG_right                    = [ betax(n+1)*SGx(n+1); betay(n+1)*SGy(n+1); betaz(n+1)*SGz(n+1) ];
    end
         
    
    A1_dense = [Wd*real(J);Wd*imag(J)];
    A2_dense = Wd1*A_seis;
    P = -[Wd*real(rn);Wd*imag(rn);Wd1*rn1;2*WG_right;sqrt(apha(n+1))*(L)*(m-mapr);sqrt(apha1(n+1))*(L)*(m1Log-m1aprLog)];
     Wm_left = blkdiag(sqrt(apha(n+1))*(L),sqrt(apha1(n+1))*(L));  
     if GramianType == 1
         if(n==1)
             lsqr_maxit = 400;
         else
             lsqr_maxit = 200;
         end
     else
             lsqr_maxit = 200;
     end
    [dm,flag,relres,iter,resvec,lsvec] = lsqr_Implicit_A_MaxIter_Ap2(A1_dense,A2_dense,Wd,Wd1,WG_left, Wm_left,P,lsqr_tol,lsqr_maxit );
    dm1  = dm(1:Nm);   
    dm2  = dm(Nm+1:end); 

    m_full   = mpre + dm1;
    sig_full = sigx;
    sig_full(minv_ind) = exp(Logm2m(m_full,mmin,mmax,pp));
    dpr  = FEM_forward(freq,ns,recpar,par,omega,np,npo,MeshTop,bnd,hl,sl,al,sig0x,sig0y,sig0z,sig_full,sig_full,...
           sig_full,Jacob,JacobInv,JacobDet,Gauss3DW,Gauss,EPGauss,compress,symFlag,fileformat,Num_Comm,Q_Ex,Q_Ey,Q_Ez,Q_Hx,Q_Hy,Q_Hz,mpirun,Xtop,Ytop,Ztop,inv_flag);
    dpre_full = dpr;
    ntest     = 100;
    kn_test   = linspace(0.3,1,ntest);
    Pa_test   = zeros(ntest,1);
    for it = 1:ntest
        dpre_test = dpr1 + kn_test(it)*(dpre_full-dpr1);
        Pa_test(it) = norm(Wd*(dpre_test-dobs));
    end
    kn = kn_test(Pa_test == min(Pa_test));
    kn = min(kn)
   % kn1  =  Step_Linear_Search_log(Wd1,dobs1,m1Log,mtimes,dm2,m1min,m1max,Nmx,Nmy,Nmz)
     kn1  = 0.2
    if GramianType == 1
        [gkn,gkn1] = Step_Linear_Search_log_joint(m,dm1,mmin,mmax,m1Log,dm2,m1min,m1max,pp)
        gkn        = 0.1;
        gkn1       = 0.1;
        kn         = (kn+gkn)/2;
        kn1        = (kn1+gkn1)/2;
        m          = m+kn*dm1;
        m1Log      = m1Log+kn1*dm2;
    else
        [gxkn,gxkn1]   = Step_Linear_Search_log_joint_G(m,dm1,mmin,mmax,m1Log,dm2,m1min,m1max,Lx,pp)
        [gykn,gykn1]   = Step_Linear_Search_log_joint_G(m,dm1,mmin,mmax,m1Log,dm2,m1min,m1max,Ly,pp)
        [gzkn,gzkn1]   = Step_Linear_Search_log_joint_G(m,dm1,mmin,mmax,m1Log,dm2,m1min,m1max,Lz,pp)
        kn             = kn/2 + gxkn/6 + gykn/6 + gzkn/6;
        kn1            = kn1/2 + gxkn1/6 + gykn1/6 + gzkn1/6;
        m              = m+kn*dm1;
        m1Log          = m1Log+kn1*dm2;
    end

    m1  =  Logm2m1(m1Log,m1min,m1max);
          
    m1  = flip(reshape(m1,Nmx,Nmy,Nmz),3);  
    m1  = m1(:);                          
    Creat_New_Vgrids(m1);
  
     if(norm(m-mpre)/norm(mpre) <= JacobUpdate_tol)
         JacobUpdate = 0;
     else
         JacobUpdate = 1;
     end
    JacobUpdate
    
    mpre = m;
        
end
