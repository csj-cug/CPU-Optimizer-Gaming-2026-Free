function [kn] = Step_Linear_Search_log(A,d,ml,dm,mmin,mmax,p)
ntest = 100;
kn_test = linspace(0.1,1.0,ntest);
Misfit_test = zeros(ntest,1);
for it = 1:ntest
    ml_test = ml+kn_test(it)*dm;
     m = Logm2m(ml_test,mmin,mmax,p);
   %  m  = ml_test;
%     m = exp(ml_test)+m_min;
    Misfit_test(it) = norm(A*m-d);
end
kn = min(kn_test(Misfit_test == min(Misfit_test)));