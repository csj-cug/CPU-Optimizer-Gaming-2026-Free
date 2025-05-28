fid = fopen('b1.txt','w');
importdata grid3dg.in;
C = ans;
format short
%   读取grid3dg里面的r、lat、long %
r1 = strsplit( cell2mat( C(17) ) ); %C17:r   cell2mat:将元胞数组转换为基础数据类型的普通数组    strsplit - 在指定分隔符处拆分字符串或字符向量
r_min =  str2num( cell2mat( r1(1) ) );        %str2num - 将字符数组或字符串转换为数值数组
r_max =  str2num( cell2mat( r1(2) ) );

lat1 = strsplit( cell2mat( C(18) ) );          %C18:lat
lat_min =  str2num( cell2mat( lat1(1) ) );
lat_max =  str2num( cell2mat( lat1(2) ) );

long1 = strsplit( cell2mat( C(19) ) );          %C19:long
long_min =  str2num( cell2mat( long1(1) ) );
long_max =  str2num( cell2mat( long1(2) ) );

r = abs(r_max-r_min);
lat = abs(lat_max-lat_min);
long = abs(long_max-long_min);


 r_scale = r/(78-1);
 lat_scale = lat/(78-1);
 long_scale = long/(24-1);


 Nmx                             = 40;  
 Nmy                             = 40;
 Nmz                             = 26;

A    = zeros(Nmx,Nmy,Nmz);


%%%%%%%%%%    设置背景速度   %%%%%%%%%%%%
  for i = 1:Nmx                         %经度   12+2  上下两个边界
      for j = 1:Nmy                     %纬度
          for k = 1:Nmz                 %深度

                      A(i,j,k) = 3.0 + 0.1*(26-k);

          end
     end
  end

 

A1 = A(:);
F = A1';
fprintf(fid,'%f\n',F);
