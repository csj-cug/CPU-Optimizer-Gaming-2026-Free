
function [kn1] = Step_Linear_Search_log(Wd1,dobs,m1Log,mtimes,dm1,m1min,m1max,Nmx,Nmy,Nmz)
ntest      = 100;
kn_test    = linspace(0.3,1,ntest);
Pa_test    = zeros(ntest,1);
dpr1       = mtimes;               %旧预测值
m0Log      = m1Log + dm1;
m0         = Logm2m1(m0Log,m1min,m1max);
m0         = flip(reshape(m0,Nmx,Nmy,Nmz),3);
m0         = m0(:);
Creat_New_Vgrids(m0);            %建立新的vgrids用于步长的正演
fm3d;                            % forward,得到走时数据arrival.dat和灵敏度矩阵frechet.dat
dpre_full  = read_mtimes;        %新预测值a
for it = 1:ntest
  %  ml_test = ml+kn_test(it)*dm;
    dpre_test = dpr1 + kn_test(it)*(dpre_full - dpr1);
 %     m = Logm2m(ml_test,mmin,mmax);
% %     m = exp(ml_test)+m_min;
    Pa_test(it) = norm(Wd1*(dpre_test-dobs));
end
kn1 = min(kn_test(Pa_test == min(Pa_test)));

         