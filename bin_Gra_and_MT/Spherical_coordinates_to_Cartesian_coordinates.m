function [dmx,dmy,dmz,Nmx,Nmy,Nmz,mx,my,mz,p,t,t_xyz_center,Rough_N_Neighbour] = Spherical_coordinates_to_Cartesian_coordinates
importdata grid3dg.in;
C = ans;
format long G 
%   读取grid3dg里面的r、lat、long %
r1    = strsplit( cell2mat( C(17) ) );       
r_min =  str2num( cell2mat( r1(1) ) );       
r_max =  str2num( cell2mat( r1(2) ) );

lat1    = strsplit( cell2mat( C(18) ) );   
lat_min =  str2num( cell2mat( lat1(1) ) );
lat_max =  str2num( cell2mat( lat1(2) ) );

long1    = strsplit( cell2mat( C(19) ) );    
long_min =  str2num( cell2mat( long1(1) ) );
long_max =  str2num( cell2mat( long1(2) ) );

r    = abs(r_max-r_min);
lat  = abs(lat_max-lat_min);
long = abs(long_max-long_min);


mesh1     =   strsplit( cell2mat( C(35) ) );     
mesh_r    =   str2num( cell2mat( mesh1(1) ) );
mesh2     =   strsplit( cell2mat( C(36) ) );      
mesh_lat  =   str2num( cell2mat( mesh2(1) ) );
mesh3     =   strsplit( cell2mat( C(37) ) );      
mesh_long =   str2num( cell2mat( mesh3(1) ) );
%网格节点数
r_node     =   mesh_r + 3;
lat_node   =   mesh_lat + 3;    
long_node  =   mesh_long + 3;

%单元个数
mesh_r     =   mesh_r + 2;
mesh_lat   =   mesh_lat + 2;
mesh_long  =   mesh_long + 2;


r_scale    =  r./mesh_r;
lat_scale  =  lat./mesh_lat;
long_scale =   long./mesh_long;

dmx = long_scale.*100.*1000;
dmx = double(vpa(dmx));
dmy = lat_scale.*100.*1000;   %小尺度不考虑曲率的情况下0.01°= 1000m
dmy = double(vpa(dmy));
dmz = r_scale*1000;            %转化单位为米
dmz = double(vpa(dmz));

lenx = dmx * mesh_long;
leny = dmy * mesh_lat;
lenz = dmz * mesh_r;

[Xp,Yp,Zp]  = meshgrid(linspace(0,lenx,long_node),linspace(0,leny,lat_node),linspace(0,lenz,r_node));
xp                  = Xp(:);
yp                  = Yp(:);
zp                  = Zp(:);
x                   = unique(xp)';
y                   = unique(yp)';
z                   = unique(zp)';
[X, Y, Z]           = ndgrid(x,y,z);
MeshReg             = connectivity_brick_ByNode(X,Y,Z);
p                   = MeshReg.p;                  %节点坐标
t                   = MeshReg.t;                  %节点编号
t_xyz_center        = MeshReg.t_xyz_center;       %中心坐标


Rough_N_Neighbour   = 26;
Nmx                 = mesh_long;
Nmy                 = mesh_lat;
Nmz                 = mesh_r;
%计算中心点坐标
mx0                             = 0;  %Starting point in x direction
my0                             = 0;
mz0                             = 0;

mx                              = mx0+dmx/2 : dmx : mx0+dmx*mesh_long-dmx/2 ;   
my                              = my0+dmy/2 : dmy : my0+dmy*mesh_lat-dmy/2 ;
mz                              = mz0+dmz/2 : dmz : mz0+dmz*mesh_r-dmz/2 ;
[my,mx,mz]                      = meshgrid(my,mx,mz);
Realm                           = [mx(:),my(:),mz(:)];  %中心点坐标





    







