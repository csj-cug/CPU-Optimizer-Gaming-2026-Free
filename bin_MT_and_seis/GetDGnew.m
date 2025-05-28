function [DG]=GetDGnew(m1,m2,WG)
     DG1                 = 2*( WG'*(WG*m1) * dot(WG*m2,WG*m2) - WG'*(WG*m2) * dot(WG*m1,WG*m2) );
     DG2                 = 2*( WG'*(WG*m2) * dot(WG*m1,WG*m1) - WG'*(WG*m1) * dot(WG*m1,WG*m2) );

    DG = [DG1;DG2]';
