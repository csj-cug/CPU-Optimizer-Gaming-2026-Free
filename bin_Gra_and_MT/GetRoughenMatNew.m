function L = GetRoughenMatNew(p,t,minv_ind,t_xyz_center,abcd_vol,N)
%Get the roughtness matrix for tethedral mesh in the inversion domain

% [idx, dist] = knnsearch(p,recpar,'dist','euclidean','k',N);
% vol = abs(abcd_vol(minv_ind,end));
vol = abs(abcd_vol);
t_xyz_center_inv = t_xyz_center(minv_ind,:);
[idx1, dist1] = knnsearch(t_xyz_center_inv,t_xyz_center_inv,'dist','euclidean','k',N+1);
idx = idx1(:,2:end);
dist=dist1(:,2:end);

ii2D = repmat((1:length(minv_ind))',[1,N+1]);
jj2D = idx1;
val2D=-ones(size(ii2D));

vol2D=vol(idx1);
vol_sum = sum(vol2D(:,2:end),2);
vol_sum2D=repmat(vol_sum,[1,N+1]);

val2D(:,2:end) = -sqrt(vol2D(:,2:end)./vol_sum2D(:,2:end)).*(1./dist1(:,2:end));
val2D(:,1) = -sum(val2D(:,2:end),2);



% val2D_tmp = 0*val2D;
% vol2D1 = repmat(vol2D(:,1),[1,N+1]);
% 
% % dist1_norm = dist1./repmat(dist1(:,end),[1,N+1]);
% % dist1_norm(:,1)=0;
% 
% % val2D_tmp(:,2:end) = sqrt(vol2D1(:,2:end)./vol_sum2D(:,2:end)).*(1./dist1_norm(:,2:end));
% val2D_tmp(:,2:end) = sqrt(vol2D1(:,2:end)./vol_sum2D(:,2:end)).*(1./dist1(:,2:end));
% val2D(:,1)=sum(val2D_tmp,2);

ii2D_T=ii2D';
jj2D_T=jj2D';
val2D_T=val2D';

L=sparse(ii2D_T(:),jj2D_T(:),val2D_T(:),length(minv_ind),length(minv_ind));