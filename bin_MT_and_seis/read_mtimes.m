function mtimes = read_mtimes
[a1,a2,a3,a4,a5,a6,a7]=textread('arrivals.dat','%f%f%f%f%f%s%s');
mtimes = a5;
