function inptdata = AddGaussNoise(inptdata,noise_lel_Z,noise_lel_T)
%% Add Gauss Noise

freq = unique(inptdata(:,5));
nfreq = length(freq);
inptdata1 = [];   
for ifreq = 1:nfreq
    %Add 2% Gauss noise to impedance data, 5% Gauss noise to tipper data
	indfZxx  = find(inptdata(:,5)==freq(ifreq) & inptdata(:,4)==3);
    indfZxy  = find(inptdata(:,5)==freq(ifreq) & inptdata(:,4)==1);
    indfZyx  = find(inptdata(:,5)==freq(ifreq) & inptdata(:,4)==2);
    indfZyy  = find(inptdata(:,5)==freq(ifreq) & inptdata(:,4)==4);
    indfTx   = find(inptdata(:,5)==freq(ifreq) & inptdata(:,4)==5);
    indfTy   = find(inptdata(:,5)==freq(ifreq) & inptdata(:,4)==6);

    Zxx      = inptdata(indfZxx,6) + 1i*inptdata(indfZxx,7); 
    Zxy      = inptdata(indfZxy,6) + 1i*inptdata(indfZxy,7);
    Zyx      = inptdata(indfZyx,6) + 1i*inptdata(indfZyx,7);
    Zyy      = inptdata(indfZyy,6) + 1i*inptdata(indfZyy,7);

%     Zxxnoise = randn(size(Zxx)).*noise_lel_Z.*sqrt(abs(Zxy.*Zyx));
%     Zxynoise = randn(size(Zxy)).*noise_lel_Z.*sqrt(abs(Zxy.*Zyx));
%     Zyxnoise = randn(size(Zyx)).*noise_lel_Z.*sqrt(abs(Zxy.*Zyx));
%     Zyynoise = randn(size(Zyy)).*noise_lel_Z.*sqrt(abs(Zxy.*Zyx));
% 
% 
%     Zxxreal  = real(Zxx) + Zxxnoise;
%     Zxximag  = imag(Zxx) + Zxxnoise;
%     Zxx      = Zxxreal + 1i*Zxximag;
% 
%     Zxyreal  = real(Zxy) + Zxynoise;
%     Zxyimag  = imag(Zxy) + Zxynoise;
%     Zxy      = Zxyreal + 1i*Zxyimag;
% 
%     Zyxreal  = real(Zyx) + Zyxnoise;
%     Zyximag  = imag(Zyx) + Zyxnoise;
%     Zyx      = Zyxreal + 1i*Zyximag;
% 
%     Zyyreal  = real(Zyy) + Zyynoise;
%     Zyyimag  = imag(Zyy) + Zyynoise;
%     Zyy      = Zyyreal + 1i*Zyyimag;

    Zxx = Zxx.*(1 + randn(size(Zxx)).*noise_lel_Z);
    Zxy = Zxy.*(1 + randn(size(Zxy)).*noise_lel_Z);
    Zyx = Zyx.*(1 + randn(size(Zyx)).*noise_lel_Z);
    Zyy = Zyy.*(1 + randn(size(Zyy)).*noise_lel_Z);

    inptdataZxx(:,1:3) = inptdata(inptdata(:,5)==freq(1) & inptdata(:,4)==3,1:3);
    inptdataZxy(:,1:3) = inptdata(inptdata(:,5)==freq(1) & inptdata(:,4)==1,1:3);
    inptdataZyx(:,1:3) = inptdata(inptdata(:,5)==freq(1) & inptdata(:,4)==2,1:3);
    inptdataZyy(:,1:3) = inptdata(inptdata(:,5)==freq(1) & inptdata(:,4)==3,1:3);
    inptdataZxx(:,4)   = 3;
    inptdataZxx(:,5)   = freq(ifreq);
    inptdataZxy(:,4)   = 1;
    inptdataZxy(:,5)   = freq(ifreq);
    inptdataZyx(:,4)   = 2;
    inptdataZyx(:,5)   = freq(ifreq);
    inptdataZyy(:,4)   = 4;
    inptdataZyy(:,5)   = freq(ifreq);

    inptdataZxx(:,6:7) = [real(Zxx),imag(Zxx)];
    inptdataZxy(:,6:7) = [real(Zxy),imag(Zxy)];
    inptdataZyx(:,6:7) = [real(Zyx),imag(Zyx)];
    inptdataZyy(:,6:7) = [real(Zyy),imag(Zyy)];

    Tx     = inptdata(indfTx,6) + 1i*inptdata(indfTx,7); 
    Ty     = inptdata(indfTy,6) + 1i*inptdata(indfTy,7);

    Tx  = Tx.*(1 + randn(size(Tx)).*noise_lel_T);
    Ty  = Ty.*(1 + randn(size(Ty)).*noise_lel_T);

    inptdataTx(:,1:3) = inptdata(inptdata(:,5)==freq(1) & inptdata(:,4)==5,1:3);
    inptdataTy(:,1:3) = inptdata(inptdata(:,5)==freq(1) & inptdata(:,4)==6,1:3);
    inptdataTx(:,4)   = 5;
    inptdataTx(:,5)   = freq(ifreq);
    inptdataTy(:,4)   = 6;
    inptdataTy(:,5)   = freq(ifreq);

    inptdataTx(:,6:7) = [real(Tx),imag(Tx)];
    inptdataTy(:,6:7) = [real(Ty),imag(Ty)];

    PTxx     = ( real(Zyy).*imag(Zxx) - real(Zxy).*imag(Zyx))./(real(Zxx).*real(Zyy) - real(Zxy).*real(Zyx));
    PTxy     = ( real(Zyy).*imag(Zxy) - real(Zxy).*imag(Zyy))./(real(Zxx).*real(Zyy) - real(Zxy).*real(Zyx));
    PTyx     = (-real(Zyx).*imag(Zxx) + real(Zxx).*imag(Zyx))./(real(Zxx).*real(Zyy) - real(Zxy).*real(Zyx));
    PTyy     = (-real(Zyx).*imag(Zxy) + real(Zxx).*imag(Zyy))./(real(Zxx).*real(Zyy) - real(Zxy).*real(Zyx));

    inptdataPTxx(:,1:3) = inptdata(inptdata(:,5)==freq(1) & inptdata(:,4)==7,1:3);
    inptdataPTyy(:,1:3) = inptdata(inptdata(:,5)==freq(1) & inptdata(:,4)==8,1:3);
    inptdataPTxy(:,1:3) = inptdata(inptdata(:,5)==freq(1) & inptdata(:,4)==9,1:3);
    inptdataPTyx(:,1:3) = inptdata(inptdata(:,5)==freq(1) & inptdata(:,4)==10,1:3);
    
    inptdataPTxx(:,4)   = 7;
    inptdataPTxx(:,5)   = freq(ifreq);
    inptdataPTyy(:,4)   = 8;
    inptdataPTyy(:,5)   = freq(ifreq);
    inptdataPTxy(:,4)   = 9;
    inptdataPTxy(:,5)   = freq(ifreq);
    inptdataPTyx(:,4)   = 10;
    inptdataPTyx(:,5)   = freq(ifreq);

    inptdataPTxx(:,6:7) = [real(PTxx),imag(PTxx)];
    inptdataPTyy(:,6:7) = [real(PTyy),imag(PTyy)];
    inptdataPTxy(:,6:7) = [real(PTxy),imag(PTxy)];
    inptdataPTyx(:,6:7) = [real(PTyx),imag(PTyx)];
    

    inptdataf = [inptdataZxy; inptdataZyx;inptdataZxx; inptdataZyy; inptdataTx; inptdataTy; inptdataPTxx; inptdataPTyy; inptdataPTxy; inptdataPTyx];    
    inptdata1 = [inptdata1;inptdataf];
end
inptdata = inptdata1;


