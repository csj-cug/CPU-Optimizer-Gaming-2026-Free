function intoutFEM = getfield_new(m,Q_Ex,Q_Ey,Q_Ez,Q_Hx,Q_Hy,Q_Hz,recpar,hl,sl,al,par,freq,omega,isrc,Xtop,Ytop,Ztop,p)

global mu0;
global sig_air;

%% Modified on 1/8,2020
%% Considering air layer
sig_air_top = abs(min(p(:,3)));
hl = [sig_air_top,hl];
sl = [sig_air,sl];
al = [1,al];

Xtop1 = Xtop(:);
Ytop1 = Ytop(:);
Ztop1 = Ztop(:);
rectp = griddata(Xtop1,Ytop1,Ztop1,recpar(:,1),recpar(:,2),'linear');
recpar(:,3) = recpar(:,3) + rectp;

recparex = recpar(recpar(:,4)==1 & recpar(:,5)==isrc & recpar(:,6)==freq,:);
recparey = recpar(recpar(:,4)==2 & recpar(:,5)==isrc & recpar(:,6)==freq,:);
recparez = recpar(recpar(:,4)==3 & recpar(:,5)==isrc & recpar(:,6)==freq,:);
recparhx = recpar(recpar(:,4)==4 & recpar(:,5)==isrc & recpar(:,6)==freq,:);
recparhy = recpar(recpar(:,4)==5 & recpar(:,5)==isrc & recpar(:,6)==freq,:);
recparhz = recpar(recpar(:,4)==6 & recpar(:,5)==isrc & recpar(:,6)==freq,:);

%get background fields
[Ep,Hp]=green3d(freq,hl,sl,al,recparex(:,1),recparex(:,2),recparex(:,3)+sig_air_top,par);

% apply MT source
if par(1) == 1 && length(par)==2
    switch par(2)
        case 1
            Ep(:,[2,3]) =0;
            Hp(:,[1,3]) =0;
        case 2
            Ep(:,[1,3]) =0;
            Hp(:,[2,3]) =0;
        otherwise
            error('par(2) must equal to 1 or 2');
    end
end

Epx = Ep(:,1);
Epy = Ep(:,2);
Epz = Ep(:,3);
Hpx = Hp(:,1);
Hpy = Hp(:,2);
Hpz = Hp(:,3);

Eax = Q_Ex*m;
Eay = Q_Ey*m;
Eaz = Q_Ez*m;

%Modified on April 10 2019
Hax = -1/(1i*omega*mu0)*(Q_Hx*m);
Hay = -1/(1i*omega*mu0)*(Q_Hy*m);
Haz = -1/(1i*omega*mu0)*(Q_Hz*m);

intoutFEMex = [recparex,real(Epx),imag(Epx),real(Eax),imag(Eax)];
intoutFEMey = [recparey,real(Epy),imag(Epy),real(Eay),imag(Eay)];
intoutFEMez = [recparez,real(Epz),imag(Epz),real(Eaz),imag(Eaz)];
intoutFEMhx = [recparhx,real(Hpx),imag(Hpx),real(Hax),imag(Hax)];
intoutFEMhy = [recparhy,real(Hpy),imag(Hpy),real(Hay),imag(Hay)];
intoutFEMhz = [recparhz,real(Hpz),imag(Hpz),real(Haz),imag(Haz)];

intoutFEM = [intoutFEMex;intoutFEMey;intoutFEMez;intoutFEMhx;intoutFEMhy;intoutFEMhz];
