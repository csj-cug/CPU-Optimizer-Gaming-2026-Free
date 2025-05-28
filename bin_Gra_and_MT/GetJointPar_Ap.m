function jointpar = GetJointPar_Ap(A1,A2,Wd1,Wd2,DG,qq,c,n)
     if n == 1
         jointpar = 0;
     else
        rng('default')
        Nm = size(A1,2);
        x = rand(2*Nm,1);
        x1 = x(1:Nm);
        x2 = x(Nm+1:end);
        norm1 = A1'*(Wd1*(Wd1*(A1*x1)));
        norm2 = A2'*(Wd2*(Wd2*(A2*x2)));
        jointpar = norm([norm1;norm2]) / norm(DG'*(DG*x));
     %     jointpar = sqrt( qq*jointpar); %Cai
        jointpar = sqrt( qq*jointpar/(n.^c) ); %Cai
    end
end