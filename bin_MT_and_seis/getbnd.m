function Bdirc   = getbnd(bnd)

%% Get boundary condition; for the current version, only consider uniform Dirchelet boundary condtion

bnd1        = [];
bnd1(:,1)   = 1:size(bnd,1);
bnd1(:,2)   = bnd;
IndBnd      = bnd1(bnd1(:,2)==1,1);       % Index of boundary nodes
numbd       = length(IndBnd);             % Number of nodes on the boundary;

Bdirc       = zeros(numbd,2);
Bdirc(:,1)  = IndBnd;
Bdirc(:,2)  = 0;