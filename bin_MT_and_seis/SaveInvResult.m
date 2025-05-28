%============================================================

InvResult.nIT             = n;
InvResult.m1              = m1;
InvResult.alpha1          = alpha1;
InvResult.misfit1         = misfit1;
InvResult.M1Inv           = M1Inv;
InvResult.Time            = toc;
InvResult.Date            = datestr(now);
InvResult.InvPar          = InvPar;
InvResult.d1pre           = mtimes;
InvResult.dobs             = dobs;
InvResult.dm              = dm;

save InvResult InvResult