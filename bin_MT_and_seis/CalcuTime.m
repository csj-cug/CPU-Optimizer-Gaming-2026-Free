function CalcuTime(toc)
hh = floor(toc/3600);
t1 = mod(toc,3600);
mm = floor(t1/60);
ss = floor(mod(t1,60));
fprintf('The code running time is: %dhours %dmins %dsecs \n',hh,mm,ss);
% fprintf('The iterations takes: %dhours %dmins %dsecs \n',hh,mm,ss);