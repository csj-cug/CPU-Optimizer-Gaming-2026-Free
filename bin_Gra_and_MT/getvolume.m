function Ve = getvolume(JacobDet,Gauss3DW)
%% Compute the volume of each element by Gauss integral
Ne      = size(JacobDet,1);
Weight  = repmat(Gauss3DW,Ne,1);
Ve      = sum(JacobDet.*Weight,2);