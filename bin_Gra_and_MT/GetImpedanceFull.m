function [inptdata,inptdataB] = GetImpedanceFull(intout)
%% Compute the impedance from the field
%% This version compute the full impedance tensor
freq    = unique(intout(:,6));
src     = unique(intout(:,5));
comp    = unique(intout(:,4));

nf      = length(freq);
ns      = length(src);
ncomp   = length(comp);


inptdata=[];%total impedance
inptdataB=[];%background impedance

for ifreq = 1:nf
    intoutf = intout(intout(:,6)==freq(ifreq),:);
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
    
    Ex1b     = (intoutf(indEx1,7)) + 1i*(intoutf(indEx1,8));
    Ey1b     = (intoutf(indEy1,7)) + 1i*(intoutf(indEy1,8));
    Hx1b     = (intoutf(indHx1,7)) + 1i*(intoutf(indHx1,8));
    Hy1b     = (intoutf(indHy1,7)) + 1i*(intoutf(indHy1,8));
    Hz1b     = (intoutf(indHz1,7)) + 1i*(intoutf(indHz1,8));
    
    Ex2b     = (intoutf(indEx2,7)) + 1i*(intoutf(indEx2,8));
    Ey2b     = (intoutf(indEy2,7)) + 1i*(intoutf(indEy2,8));
    Hx2b     = (intoutf(indHx2,7)) + 1i*(intoutf(indHx2,8));
    Hy2b     = (intoutf(indHy2,7)) + 1i*(intoutf(indHy2,8));
    Hz2b     = (intoutf(indHz2,7)) + 1i*(intoutf(indHz2,8));
    
    Zxx      = (Ex1.*Hy2 - Ex2.*Hy1)./(Hx1.*Hy2-Hx2.*Hy1);
    Zxxb     = (Ex1b.*Hy2b - Ex2b.*Hy1b)./(Hx1b.*Hy2b-Hx2b.*Hy1b);
    
    Zxy      = (Ex2.*Hx1 - Ex1.*Hx2)./(Hx1.*Hy2 - Hx2.*Hy1);
    Zxyb     = (Ex2b.*Hx1b - Ex1b.*Hx2b)./(Hx1b.*Hy2b - Hx2b.*Hy1b);
    
    Zyx      = (Ey1.*Hy2 - Ey2.*Hy1)./(Hx1.*Hy2-Hx2.*Hy1);
    Zyxb     = (Ey1b.*Hy2b - Ey2b.*Hy1b)./(Hx1b.*Hy2b-Hx2b.*Hy1b);
    
    Zyy      = (Ey2.*Hx1 - Ey1.*Hx2)./(Hx1.*Hy2-Hx2.*Hy1);
    Zyyb     = (Ey2b.*Hx1b - Ey1b.*Hx2b)./(Hx1b.*Hy2b-Hx2b.*Hy1b);
    
    %Tipper
    Tx       = (Hz1.*Hy2 - Hz2.*Hy1)./(Hx1.*Hy2-Hx2.*Hy1);
    Txb      = (Hz1b.*Hy2b - Hz2b.*Hy1b)./(Hx1b.*Hy2b-Hx2b.*Hy1b);
    Ty       = (Hz2.*Hx1 - Hz1.*Hx2)./(Hx1.*Hy2-Hx2.*Hy1);
    Tyb      = (Hz2b.*Hx1b - Hz1b.*Hx2b)./(Hx1b.*Hy2b-Hx2b.*Hy1b);
    
    %Phase Tensor
    %refered to Patro, 2013
    
    PTxx     = ( real(Zyy).*imag(Zxx) - real(Zxy).*imag(Zyx))./(real(Zxx).*real(Zyy) - real(Zxy).*real(Zyx));
    PTxxb    = ( real(Zyyb).*imag(Zxxb) - real(Zxyb).*imag(Zyxb))./(real(Zxxb).*real(Zyyb) - real(Zxyb).*real(Zyxb));
    
    PTxy     = ( real(Zyy).*imag(Zxy) - real(Zxy).*imag(Zyy))./(real(Zxx).*real(Zyy) - real(Zxy).*real(Zyx));
    PTxyb    = ( real(Zyyb).*imag(Zxyb) - real(Zxyb).*imag(Zyyb))./(real(Zxxb).*real(Zyyb) - real(Zxyb).*real(Zyxb));
    
    PTyx     = (-real(Zyx).*imag(Zxx) + real(Zxx).*imag(Zyx))./(real(Zxx).*real(Zyy) - real(Zxy).*real(Zyx));
    PTyxb    = (-real(Zyxb).*imag(Zxxb) + real(Zxxb).*imag(Zyxb))./(real(Zxxb).*real(Zyyb) - real(Zxyb).*real(Zyxb));
    
    PTyy     = (-real(Zyx).*imag(Zxy) + real(Zxx).*imag(Zyy))./(real(Zxx).*real(Zyy) - real(Zxy).*real(Zyx));
    PTyyb    = (-real(Zyxb).*imag(Zxyb) + real(Zxxb).*imag(Zyyb))./(real(Zxxb).*real(Zyyb) - real(Zxyb).*real(Zyxb));
    
    
    %% Compute Zxy
    inptdataZxy = zeros(size(Ex1,1),7);
    inptdataZxy(:,1:3)=intoutf(indEx1,1:3);
    inptdataZxy(:,4)=1;
    inptdataZxy(:,5)=freq(ifreq);
    inptdataZxy(:,6:7)=[real(Zxy),imag(Zxy)];
    
    inptdataZxyB=inptdataZxy;
    inptdataZxyB(:,6:7)=[real(Zxyb),imag(Zxyb)];
    
    
    %% Compute Zyx
    inptdataZyx=inptdataZxy;
    inptdataZyx(:,4)=2;
    inptdataZyx(:,6:7)=[real(Zyx),imag(Zyx)];
    
    inptdataZyxB=inptdataZyx;
    inptdataZyxB(:,6:7)=[real(Zyxb),imag(Zyxb)];
    
    %% Compute Zxx
    inptdataZxx = inptdataZxy;
    inptdataZxx(:,4)=3;
    inptdataZxx(:,5)=freq(ifreq);
    inptdataZxx(:,6:7)=[real(Zxx),imag(Zxx)];
    
    inptdataZxxB=inptdataZxx;
    inptdataZxxB(:,6:7)=[real(Zxxb),imag(Zxxb)];
    
    %% Compute Zyy
    inptdataZyy=inptdataZxy;
    inptdataZyy(:,4)=4;
    inptdataZyy(:,6:7)=[real(Zyy),imag(Zyy)];
    
    inptdataZyyB=inptdataZyy;
    inptdataZyyB(:,6:7)=[real(Zyyb),imag(Zyyb)];
    
    %% Compute Tx
    inptdataTx = zeros(size(Ex1,1),7);
    inptdataTx(:,1:3)=intoutf(indEx1,1:3);
    inptdataTx(:,4)=5;
    inptdataTx(:,5)=freq(ifreq);
    inptdataTx(:,6:7)=[real(Tx),imag(Tx)];
    
    inptdataTxB=inptdataTx;
    inptdataTxB(:,6:7)=[real(Txb),imag(Txb)];
    
    %% Compute Ty
    inptdataTy = inptdataTx;
    inptdataTy(:,4)=6;
    inptdataTy(:,6:7)=[real(Ty),imag(Ty)];
    
    inptdataTyB=inptdataTy;
    inptdataTyB(:,6:7)=[real(Tyb),imag(Tyb)];
    
    %% Compute PTxx
    inptdataPTxx = zeros(size(Ex1,1),6);
    inptdataPTxx(:,1:3)=intoutf(indEx1,1:3);
    inptdataPTxx(:,4)=7;
    inptdataPTxx(:,5)=freq(ifreq);
    inptdataPTxx(:,6:7)=[real(PTxx),imag(PTxx)];
    
    inptdataPTxxB=inptdataPTxx;
    inptdataPTxxB(:,6:7)=[real(PTxxb),imag(PTxxb)];
    
    %% Compute PTyy
    inptdataPTyy = inptdataPTxx;
    inptdataPTyy(:,4)=8;
    inptdataPTyy(:,6:7)=[real(PTyy),imag(PTyy)];
    
    inptdataPTyyB=inptdataPTyy;
    inptdataPTyyB(:,6:7)=[real(PTyyb),imag(PTyyb)];
    
    %% Compute PTxy
    inptdataPTxy = inptdataPTxx;
    inptdataPTxy(:,4)=9;
    inptdataPTxy(:,6:7)=[real(PTxy),imag(PTxy)];
    
    inptdataPTxyB=inptdataPTxy;
    inptdataPTxyB(:,6:7)=[real(PTxyb),imag(PTxyb)];
    
    %% Compute PTyx
    inptdataPTyx = inptdataPTxx;
    inptdataPTyx(:,4)=10;
    inptdataPTyx(:,6:7)=[real(PTyx),imag(PTyx)];
    
    inptdataPTyxB=inptdataPTyx;
    inptdataPTyxB(:,6:7)=[real(PTyxb),imag(PTyxb)];
    
    
    %% write data to file
    %inptdataf = [inptdataZxx; inptdataZxy;inptdataZyx; inptdataZyy];
    %inptdataBf = [inptdataZxxB; inptdataZxyB; inptdataZyxB; inptdataZyyB];

    inptdataf = [inptdataZxy; inptdataZyx;inptdataZxx; inptdataZyy; inptdataTx; inptdataTy; inptdataPTxx; inptdataPTyy; inptdataPTxy; inptdataPTyx];
    inptdataBf = [inptdataZxyB; inptdataZyxB; inptdataZxxB; inptdataZyyB; inptdataTxB; inptdataTyB; inptdataPTxxB; inptdataPTyyB; inptdataPTxyB; inptdataPTyxB];
    
    inptdata = [inptdata;inptdataf];
    inptdataB = [inptdataB;inptdataBf];
 
end

1;
