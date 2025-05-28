 function A_seis = frechet(Nd,Nm,Nmx,Nmy,Nmz)

A = dlmread('frechet.dat');
n = find(A(:,5)~=0);                
 
K = zeros(Nd,Nm);
for i=1:Nd-1
K(i,A(n(i)+1:n(i+1)-1,1))=A(n(i)+1:n(i+1)-1,2);
end
K(Nd,A(n(Nd)+1:end,1))=A(n(Nd)+1:end,2);
A_seis = reshape(K,Nd,Nmz,Nmy,Nmx);
A_seis = permute(A_seis, [1 4 3 2]);
A_seis  = flip(A_seis,4);         %new
A_seis = A_seis(:);
A_seis = reshape(A_seis,Nd,Nmx*Nmy*Nmz);




