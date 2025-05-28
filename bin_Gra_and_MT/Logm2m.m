function m = Logm2m(mLog,mmin,mmax,p)
    %m              = mmax-(mmax-mmin)./(1+exp(mLog));
 %m              = (mmin+mmax*exp(mLog))./(1+exp(mLog));
m              = (mmin+mmax*exp(p*mLog))./(1+exp(p*mLog));