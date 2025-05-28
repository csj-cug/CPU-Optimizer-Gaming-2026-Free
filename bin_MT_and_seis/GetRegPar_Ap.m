% Calculate regularization parameter using spectral method  使用谱法计算正则化参数
function  regpar  = GetRegPar_Ap(A,p,Wd,L,qq,c,n)
    rng('default')
    Nm = size(A,2);            
    x = rand(Nm,1);
    regpar = norm(p.*(A'*(Wd*(Wd*(A*(p.*x)))))) / norm(L'*(L*x));
    regpar = sqrt( qq*regpar/(n.^c) );

    