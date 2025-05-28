function [Xtop,Ytop,Ztop] = XYZTop(x,y)


%
nx=length(x);
ny=length(y);

[Xtop,Ytop]=ndgrid(x,y);

Ztop=zeros(nx,ny);



