function  [vgrids_ref,vgrids_true] = read_VgridApr 
A = dlmread('vgridstrue.in');

B = dlmread('vgridsref.in');

r_node    =   A(2,1);
lat_node  =   A(2,2);
long_node =   A(2,3);

for i = 1: r_node*lat_node*long_node
    vgrids_true(i) = A(i+4,1);
    vgrids_ref(i)  = B(i+4,1);
end
vgrids_true = vgrids_true';
vgrids_ref = vgrids_ref';

save vgrids_true vgrids_true
