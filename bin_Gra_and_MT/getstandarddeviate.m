function [inptdata1,inptdataB1] = getstandarddeviate(inptdata,inptdataB,errZoff,errZon)
%% compute measurement error or standard deviation each sites
%% 2021.04.20

freq  = unique(inptdata(:,5));
nfreq = length(freq);

inptdata1 = [];
for ifreq = 1:nfreq
    inptf   = inptdata(inptdata(:,5) == freq(ifreq),:);
    inptZxy = inptf(inptf(:,4) == 1,:);
    inptZyx = inptf(inptf(:,4) == 2,:);
    inptZxx = inptf(inptf(:,4) == 3,:);
    inptZyy = inptf(inptf(:,4) == 4,:);
    inptTx  = inptf(inptf(:,4) == 5,:);
    inptTy  = inptf(inptf(:,4) == 6,:);
    inptPTxy= inptf(inptf(:,4) == 7,:);
    inptPTyx= inptf(inptf(:,4) == 8,:);
    inptPTxx= inptf(inptf(:,4) == 9,:);
    inptPTyy= inptf(inptf(:,4) == 10,:);
    Zxytp   = inptZxy(:,end-1) + 1i*inptZxy(:,end);
    Zyxtp   = inptZyx(:,end-1) + 1i*inptZyx(:,end);
    Txtp    = inptTx(:,end-1) + 1i*inptTx(:,end);
    Tytp    = inptTy(:,end-1) + 1i*inptTy(:,end);
    PTxxtp  = inptPTxx(:,end-1) + 1i*inptPTxx(:,end);
    PTyytp  = inptPTyy(:,end-1) + 1i*inptPTyy(:,end);
    errZ    = sqrt(abs(Zxytp.*Zyxtp)); %
    errT    = 0.03 + zeros(length(Txtp),1); %refered to Kordy, 2016
    errPT   = sqrt(abs(PTxxtp.*PTyytp));
    
    %% Zxy
    inptdataZxy = zeros(size(inptZxy,1),8);
    inptdataZxy(:,1:7) = inptZxy;
    inptdataZxy(:,8)   = errZoff*errZ;
    
    %% Zyx
    inptdataZyx = zeros(size(inptZyx,1),8);
    inptdataZyx(:,1:7) = inptZyx;
    inptdataZyx(:,8)   = errZoff*errZ;
    
    %% Zxx
    inptdataZxx = zeros(size(inptZxx,1),8);
    inptdataZxx(:,1:7) = inptZxx;
    inptdataZxx(:,8)   = errZon*errZ;
    
    %% Zyy
    inptdataZyy = zeros(size(inptZyy,1),8);
    inptdataZyy(:,1:7) = inptZyy;
    inptdataZyy(:,8)   = errZon*errZ;
    
    %% Tx
    inptdataTx = zeros(size(inptTx,1),8);
    inptdataTx(:,1:7) = inptTx;
    inptdataTx(:,8)   = errT;
    
    %% Ty
    inptdataTy = zeros(size(inptTy,1),8);
    inptdataTy(:,1:7) = inptTy;
    inptdataTy(:,8)   = errT;
    
    %% PTxy
    inptdataPTxy = zeros(size(inptPTxy,1),8);
    inptdataPTxy(:,1:7) = inptPTxy;
    inptdataPTxy(:,8)   = errPT;
    
    %% PTyx
    inptdataPTyx = zeros(size(inptPTyx,1),8);
    inptdataPTyx(:,1:7) = inptPTyx;
    inptdataPTyx(:,8)   = errPT;
    
    %% PTxx
    inptdataPTxx = zeros(size(inptPTxx,1),8);
    inptdataPTxx(:,1:7) = inptPTxx;
    inptdataPTxx(:,8)   = errPT;
    
    %% PTyy
    inptdataPTyy = zeros(size(inptPTyy,1),8);
    inptdataPTyy(:,1:7) = inptPTyy;
    inptdataPTyy(:,8)   = errPT;
    
    inptdataf = [inptdataZxy; inptdataZyx;inptdataZxx; inptdataZyy; inptdataTx; inptdataTy; inptdataPTxy; inptdataPTyx; inptdataPTxx; inptdataPTyy];
    inptdata1 = [inptdata1;inptdataf];
   
end
inptdataB1 = [inptdataB,ones(size(inptdataB,1),1)];
