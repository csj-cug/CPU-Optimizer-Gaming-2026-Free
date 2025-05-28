function Wd             = GetWd1(d1,Nd)
        dobs = abs(d1);
        r    = 0.02;                %noise
        eps  = 0.0001;
% ======================================================
        Wd1 = 1./(dobs*r+eps);
        
        Wd = spdiags(Wd1,0,Nd,Nd);