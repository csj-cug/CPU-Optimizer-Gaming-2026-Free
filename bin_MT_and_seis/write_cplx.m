 function write_cplx(filename,data)

 % Writes complex dataset to ascii file

 fid = fopen(filename,'w');
 fprintf(fid,'%10.10g %10.10g\n',[real(data)';imag(data)']);
 fclose(fid);
