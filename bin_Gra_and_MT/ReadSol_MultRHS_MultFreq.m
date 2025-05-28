function m = ReadSol_MultRHS_MultFreq(ns,filename,fileformat,compress)
%Default name for output file is 'Output.txt'
if (fileformat == 0)
 fid=fopen(filename);
 [A,count] = fscanf(fid,' (%e,%e)',[2,inf]);
 A = A';
 x = A(:,1) + 1i*A(:,2);
 fclose(fid);
else
 if (~H5Z.filter_avail('H5Z_FILTER_DEFLATE'))
    disp('HDF5 compression not available, compressed file read will fail');
 end
 fid = H5F.open(filename);
 x = h5read_cplx(fid,'x_comm',1);
 H5F.close(fid);
end
 
nm = length(x)/ns;
m  = zeros(nm,ns);

for isrc = 1:ns
    m(:,isrc) = x((isrc-1)*nm+1:isrc*nm);
end