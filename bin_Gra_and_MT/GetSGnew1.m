function SG = GetSGnew(m1,m2,WG)

    mg1             = WG*m1;
    mg2             = WG*m2;
    SG              = dot(mg1,mg1)*dot(mg2,mg2)-dot(mg1,mg2)^2;