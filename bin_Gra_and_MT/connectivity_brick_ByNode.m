function output = connectivity_brick_ByNode(X,Y,Z)
%% Number scheme for the edge: First, all edges along x+ direction will be numbered;
%% Second, all edges along y+ direction will be numbered
%% Finally, all edges along z+ direction will be numbered


nx = size(X,1);   %27
ny = size(X,2);   %27
nz = size(X,3);   %40
Ne = (nx-1)*(ny-1)*(nz-1);% number of element  26364 网格数量
xn = X(:);
yn = Y(:);
zn = Z(:);


p = [xn,yn,zn];  %节点坐标
t = zeros(Ne,8);  %每个节点的编号
Num = reshape(1:size(p,1),nx,ny,nz); %Similary as X,Y,Z, but store the corresponding node number for X,Y,Z


%% Get the first node for each bricks which corresponding min(x),min(y),min(z) in each brick
NumTmp = Num(1:end-1,1:end-1,1:end-1);
t(:,1) = NumTmp(:);

%% Get the second node for each bricks which corresponding max(x),min(y),min(z) in each brick
NumTmp = Num(2:end,1:end-1,1:end-1);
t(:,2) = NumTmp(:);


%% Get the third node for each bricks which corresponding max(x),max(y),min(z) in each brick
NumTmp = Num(2:end,2:end,1:end-1);
t(:,3) = NumTmp(:);


%% Get the fourth node for each bricks which corresponding min(x),max(y),min(z) in each brick
NumTmp = Num(1:end-1,2:end,1:end-1);
t(:,4) = NumTmp(:);


%% Get the fifth node for each bricks which corresponding min(x),min(y),max(z) in each brick
NumTmp = Num(1:end-1,1:end-1,2:end);
t(:,5) = NumTmp(:);

%% Get the sixth node for each bricks which corresponding max(x),min(y),max(z) in each brick
NumTmp = Num(2:end,1:end-1,2:end);
t(:,6) = NumTmp(:);


%% Get the seventh node for each bricks which corresponding max(x),max(y),max(z) in each brick
NumTmp = Num(2:end,2:end,2:end);
t(:,7) = NumTmp(:);


%% Get the eight node for each bricks which corresponding min(x),max(y),max(z) in each brick
NumTmp = Num(1:end-1,2:end,2:end);
t(:,8) = NumTmp(:);





%% The total number of edges is the sum of edges in x direction, y direction and z direction.
%% the edge direction is defined as positive x y or z;
nedx = (nx-1)*ny*nz;
nedy = nx*(ny-1)*nz;
nedz = nx*ny*(nz-1);
ned  = nedx + nedy + nedz; %total number of edges 边总数
NumEdx1 = Num(1:end-1,:,:); %number of the first node of each x oriented edge
NumEdx2 = Num(2:end,:,:); %number of the second node of each x oriented edge
edx     = [NumEdx1(:),NumEdx2(:)];

NumEdy1 = Num(:,1:end-1,:); %number of the first node of each x oriented edge
NumEdy2 = Num(:,2:end,:); %number of the second node of each x oriented edge
edy     = [NumEdy1(:),NumEdy2(:)];


NumEdz1 = Num(:,:,1:end-1); %number of the first node of each x oriented edge
NumEdz2 = Num(:,:,2:end); %number of the second node of each x oriented edge
edz     = [NumEdz1(:),NumEdz2(:)];

ed_tmp      = [edx;edy;edz];
[ed Ind]     = sortrows(ed_tmp,[1 2]);   %对矩阵行或表行进行排序

ed_tmp_comp = ed_tmp(:,1)+1i*ed_tmp(:,2);
ed_comp     = ed(:,1)+1i*ed(:,2);

[tmp,~,Ind1] = intersect(ed_tmp_comp,ed_comp,'stable');  %%Here ed_tmp=ed(Ind1);  %polyshape 对象的交集

%% Get TE: Brick-edge connectivity
TE      = zeros(Ne,12);         %12条边编号

%Get the first edge of TE
NumEdx      = reshape(1:nedx, nx-1,ny,nz);
NumEdx_E1   = NumEdx(:,1:end-1,1:end-1);
TE(:,1)     = Ind1(NumEdx_E1(:));          

%Get the first edge of TE
NumEdx_E2   = NumEdx(:,2:end,1:end-1);
TE(:,2)     = Ind1(NumEdx_E2(:));

%Get the third edge of TE
NumEdx_E3   = NumEdx(:,1:end-1,2:end);
TE(:,3)     = Ind1(NumEdx_E3(:));

%Get the fourth edge of TE
NumEdx_E4   = NumEdx(:,2:end,2:end);
TE(:,4)     = Ind1(NumEdx_E4(:));



%Get the fifth edge of TE
NumEdy      = reshape(nedx+1:nedx+nedy, nx,ny-1,nz);
NumEdy_E5   = NumEdy(1:end-1,:,1:end-1);
TE(:,5)     = Ind1(NumEdy_E5(:));


%Get the sixth edge of TE
NumEdy_E6   = NumEdy(1:end-1,:,2:end);
TE(:,6)     = Ind1(NumEdy_E6(:));

%Get the seventh edge of TE
NumEdy_E7   = NumEdy(2:end,:,1:end-1);
TE(:,7)     = Ind1(NumEdy_E7(:));

%Get the eight edge of TE
NumEdy_E8   = NumEdy(2:end,:,2:end);
TE(:,8)     = Ind1(NumEdy_E8(:));

%Get the ninth edge of TE
NumEdz      = reshape(nedx+nedy+1:nedx+nedy+nedz, nx,ny,nz-1);
NumEdz_E9   = NumEdz(1:end-1,1:end-1,:);
TE(:,9)     = Ind1(NumEdz_E9(:));


%Get the tenth edge of TE
NumEdz_E10   = NumEdz(2:end,1:end-1,:);
TE(:,10)     = Ind1(NumEdz_E10(:));

%Get the eventhth edge of TE
NumEdz_E11   = NumEdz(1:end-1,2:end,:);
TE(:,11)     = Ind1(NumEdz_E11(:));

%Get the twelveth edge of TE
NumEdz_E12   = NumEdz(2:end,2:end,:);
TE(:,12)     = Ind1(NumEdz_E12(:));


% TE1 = GetTE(t,ed);

% dif  = unique(TE-TE1);
% dif  = dif(:);
% 
% if length(dif) ~=1
%     error('TE matrix is different for two method');
% end
% 
% if length(dif) ==1 && dif ~=0
%     error('TE matrix is different for two method');
% end


orientx = unique((xn(edx(:,2))>xn(edx(:,1))));
if (length(orientx) ~=1) || (length(orientx) ==1 && orientx ~=1)
    error('Edge along x direction should point from small x to big x');
end

orienty = unique((yn(edy(:,2))>yn(edy(:,1))));
if (length(orienty) ~=1) || (length(orienty) ==1 && orienty ~=1)
    error('Edge along y direction should point from small y to big y');
end




xyz1            = p(ed(:,1),:); %coordinate for the first nodes of each edge 每条边的第一个节点的坐标
xyz2            = p(ed(:,2),:);
ed_center       = (xyz1+xyz2)/2;
xyzdif          = xyz2-xyz1;
ed_length       = sqrt(xyzdif(:,1).^2 + xyzdif(:,2).^2 + xyzdif(:,3).^2);

ed_unit = xyzdif./[ed_length ,ed_length ,ed_length ];




xe = zeros(Ne,8);
ye = zeros(Ne,8);
ze = zeros(Ne,8);
t_xyz_center = zeros(Ne,3);
t_xyz_mm=zeros(Ne,6);
for it = 1:8
    xe(:,it) = p(t(:,it),1);
    ye(:,it) = p(t(:,it),2);
    ze(:,it) = p(t(:,it),3);
end


t_xyz_center(:,1) = sum(xe,2)/8;
t_xyz_center(:,2) = sum(ye,2)/8;
t_xyz_center(:,3) = sum(ze,2)/8;

t_xyz_mm(:,1)=min(xe,[],2);          %矢量有限元 边头尾两点坐标
t_xyz_mm(:,2)=max(xe,[],2);
t_xyz_mm(:,3)=min(ye,[],2);
t_xyz_mm(:,4)=max(ye,[],2);
t_xyz_mm(:,5)=min(ze,[],2);
t_xyz_mm(:,6)=max(ze,[],2);


TE_length = zeros(Ne,12);         

for it = 1:12
    TE_length(:,it) = ed_length(TE(:,it));
end

1;


output.p            = p;
output.t            = t;
output.t_xyz_center = t_xyz_center;
output.t_xyz_mm     = t_xyz_mm;
output.ed           = ed;
output.ed_center    = ed_center;
output.ed_length    = ed_length;
output.ed_unit      = ed_unit;
output.edx          = edx;
output.edy          = edy;
output.edz          = edz;
output.nedx         = nedx;
output.nedy         = nedy;
output.nedz         = nedz;
output.TE           = TE;
output.TE_length    = TE_length;
output.Ind1         = Ind1;

