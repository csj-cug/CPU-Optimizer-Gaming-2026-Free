function [Gauss, GaussW, Gauss3Dx, Gauss3Dy, Gauss3Dz, Gauss3DW] = GetGaussian(order)
%% function to get gaussian point and weight.
%Gauss: Gauss point
%GaussW: Gauss weight


if order == 2
    Gauss  = [-1/sqrt(3) 1/sqrt(3)];
    GaussW = [1 1];
elseif order == 3
    Gauss  = [-sqrt(3/5) 0 sqrt(3/5)];
    GaussW = [5/9 8/9 5/9];
elseif order == 4
    Gauss  = [-sqrt((3+2*sqrt(6/5))/7) -sqrt((3-2*sqrt(6/5))/7) sqrt((3-2*sqrt(6/5))/7) sqrt((3+2*sqrt(6/5))/7) ];
    GaussW = [(18-sqrt(30))/36 (18+sqrt(30))/36 (18+sqrt(30))/36 (18-sqrt(30))/36];
else
    error('No such option');
end

NGauss    = length(Gauss);
Gauss3Dx  = zeros(1,NGauss^3);
Gauss3Dy  = zeros(1,NGauss^3);
Gauss3Dz  = zeros(1,NGauss^3);
Gauss3DW  = zeros(1,NGauss^3);

for iz = 1:NGauss
    for iy = 1:NGauss
        for ix = 1:NGauss
            ind = (iz-1)*NGauss^2 + (iy-1)*NGauss+ix;
            Gauss3Dx(ind) = Gauss(ix);
            Gauss3Dy(ind) = Gauss(iy);
            Gauss3Dz(ind) = Gauss(iz);
            Gauss3DW(ind) = GaussW(ix)*GaussW(iy)*GaussW(iz);
        end
    end
end

