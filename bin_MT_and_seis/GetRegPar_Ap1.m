% Calculate regularization parameter using spectral method  使用谱法计算正则化参数
function  regpar  = GetRegPar_Ap(A,Wd,L,qq,c,n)
% A = A_seis;
% p  = p1;
% Wd = Wd1;
% L  = Wm1;
% qq = InvPar.qalpha1;
    rng('default')
    Nm = size(A,2);              %80×80×26
    x = rand(Nm,1);
    regpar = norm(A'*(Wd*(Wd*(A*x)))) / norm(L'*(L*x));
    regpar = ( qq*regpar/(n.^c) );
 %regpar = sqrt( qq*regpar/(n) );
 %regpar = sqrt( qq*regpar );

    