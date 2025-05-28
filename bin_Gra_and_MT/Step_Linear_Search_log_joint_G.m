function [kn1,kn2] = Step_Linear_Search_log_joint_G(ml1,dm1,m1_min,m1_max,ml2,dm2,m2_min,m2_max,WG,pp)
ntest = 50;
kn_test = linspace(0.1,1,ntest);
Misfit_test = zeros(ntest,ntest);
for it = 1:ntest
    ml1_test = ml1+kn_test(it)*dm1;
     m1 = exp(Logm2m(ml1_test,m1_min,m1_max,pp));
     for it1 = 1:ntest
        ml2_test = ml2+kn_test(it1)*dm2;
        m2 = Logm2m1(ml2_test,m2_min,m2_max);
        Misfit_test(it,it1) = GetSGnew1(m1,m2,WG);
     end
end
[n1,n2] = find(Misfit_test == min(Misfit_test(:)));
kn1 = kn_test(n1);
kn2 = kn_test(n2);

