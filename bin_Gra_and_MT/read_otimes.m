function otimes = read_otimes
A11 = dlmread('otimes.dat');
A111 = A11(1,1);
for i = 1:A111
    otimes(i) = A11(2*i,5);
end
otimes = otimes';