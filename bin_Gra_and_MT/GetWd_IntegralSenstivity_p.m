function Wd             = GetWd_IntegralSenstivity(J)

% ======================================================
        Nd              = size(J,1);
        Wd_diag         = zeros(Nd,1);
        for nd=1:Nd
            Wd_diag(nd)     = ((J(nd,:))*J(nd,:)')^(-0.25);
        end
        
        Wd              = spdiags(Wd_diag,0,Nd,Nd);
% ======================================================