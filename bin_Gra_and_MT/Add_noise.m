function otimes = Add_noise(noise)
[a1,a2,a3,a4,a5,a6,a7]=textread('arrivals.dat','%f%f%f%f%f%s%s');
noise  = 0.02;  
mtimes = a5;
otimes = mtimes.*(1 + randn(size(mtimes)).*noise);
save otimes otimes




