function t_xyz_center                     = Get_txyz_center(p,t)

Ne      = size(t,1);
t_xyz_center    = zeros(Ne,3);
t_xyz_center(:,1) = (p(t(:,1),1) + p(t(:,2),1) + p(t(:,3),1) + p(t(:,4),1))/4;
t_xyz_center(:,2) = (p(t(:,1),2) + p(t(:,2),2) + p(t(:,3),2) + p(t(:,4),2))/4;
t_xyz_center(:,3) = (p(t(:,1),3) + p(t(:,2),3) + p(t(:,3),3) + p(t(:,4),3))/4;