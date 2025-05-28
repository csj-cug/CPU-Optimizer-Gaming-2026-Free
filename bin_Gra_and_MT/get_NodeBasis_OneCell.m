function L_cell=get_NodeBasis_OneCell(Gauss3Dx,Gauss3Dy,Gauss3Dz)
%% get the value of node basis function in one cell at the Gauss point
%% define node xyz coordinate for the reference cell
chi_node  = [-1 1 1 -1 -1 1 1 -1];
eta_node  = [-1 -1 1 1 -1 -1 1 1];
zeta_node = [-1 -1 -1 -1 1 1 1 1];

nGauss = length(Gauss3Dx);%number of Gauss point
L_cell=zeros(nGauss,8);

for iGauss = 1:nGauss
    L_cell(iGauss,:) = (1+chi_node*Gauss3Dx(iGauss)).*(1+eta_node*Gauss3Dy(iGauss)).*(1+zeta_node*Gauss3Dz(iGauss))/8;
end