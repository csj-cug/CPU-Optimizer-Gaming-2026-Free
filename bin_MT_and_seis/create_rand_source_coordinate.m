fid = fopen('sourcesle.in','w');
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
tolerance = 0.05;

%公式 r = a + (b-a)*rand(m,n)
%其中[a,b]是范围，[m,n]是生成的数据形状。比如我想生成[-5,5]之内10个随机数
%a = -5;
%b = 5;
%r = a + (b-a) * rand(10,1);
source = 50;                  %震源个数
source_lat  = lat_min + (lat_max-lat_min) * rand(source,1);
source_long = long_min + (long_max-long_min) * rand(source,1);
source_r    = abs(r_min + (r_max-r_min) * rand(source,1));             %%注意这里r应为正数

A1 = abs(r_max)-source_r;
for i = 1:source
if (lat_min-source_lat(i) <= tolerance)
    source_lat(i) = source_lat(i)-tolerance;
end
end
for i = 1:source
if (source_lat(i)-lat_max <= tolerance)
    source_lat(i) = source_lat(i)+tolerance;
end
end
for i = 1:source
 if(source_long(i)-long_min <= tolerance)
    source_long(i) = source_long(i)+tolerance;
 end
end
for i = 1:source
 if (long_max-source_long(i) <= tolerance)
    source_long(i) = source_long(i)-tolerance;
 end
end



for i = 1:source
 if(source_r(i)-r_min <= tolerance)
    source_r(i) = source_r(i)+tolerance;
 end
end
for i = 1:source
 if (abs(r_max)-source_r(i) <= tolerance)
    source_r(i) = source_r(i)-tolerance;
 end
end
A2 = abs(r_max)-source_r;
source1 = [source_lat source_long source_r];
source1 = source1';
fprintf(fid,'%d\n',source);
fprintf(fid,'%f\t%f\t%f\n',source1);


 
