function ShowInv_Casez3
%loadInvResult
load  D:\Gauss_Newton_Inversion_for_Seismic\gauss_newton_inversion40×40×26\CaseK1\SepInv4\InvResult.mat
LoadInvResult
Sx                         = [];
Sy                         = [20];                 %切片
Sz                         = [5];             

mx = unique(mx);
my = unique(my);
mz = unique(mz);
mx = mx - 20000;
my = my - 20000;
mz = mz - 20000;
Nmx = 40;
Nmy = 40;
Nmz = 26;
ShowInvModelSlice(Sx,Sy,Sz,...
                    m1,Nmx,Nmy,Nmz,dmx,dmy,dmz,mx,my,mz)
                
 %ShowMisfit





