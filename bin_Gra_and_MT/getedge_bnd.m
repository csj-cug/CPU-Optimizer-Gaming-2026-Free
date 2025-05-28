function bnd    = getedge_bnd(ed_center)
x = ed_center(:,1);
y = ed_center(:,2);
z = ed_center(:,3);
bnd = (abs(x-min(x)) < 1) | (abs(x-max(x)) < 1) | (abs(y-min(y)) < 1) | (abs(y-max(y)) < 1) | (abs(z-min(z)) < 1) | (abs(z-max(z)) < 1);