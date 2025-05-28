function [Lx,Ly,Lz]=MLSI3D_Rec_N_InvAnalytical(p,recpar,N)

%% http://www.mathworks.com/help/curvefit/least-squares-fitting.html

% p: the node information for all nodes.
% ind: index of node where derivative need to be computed by MLSI3D
% f:    the value of physical field on the node.
% N:    nearest N number of points

Ndi = size(recpar,1);
Nm  = size(p,1);
% ABCD   = zeros(Ndi,4);
% B   = zeros(Ndi,1);
% C   = zeros(Ndi,1);
% D   = zeros(Ndi,1);

% Lx = sparse(Ndi,Nm);

ii2D = repmat((1:Ndi),N,1);
ii = ii2D(:);
jj = zeros(size(ii));
valx=zeros(size(ii));
valy=zeros(size(ii));
valz=zeros(size(ii));

Rhs = zeros(4,N,Ndi);
A1  = zeros(4,4,Ndi);
% tic;
%Most time consuming part
[idx, dist] = knnsearch(p,recpar,'dist','euclidean','k',N);
X     = ones(N,4);
for it = 1:Ndi
%     rd = sqrt((p(:,1)-recpar(it,1)).^2 + (p(:,2)-recpar(it,2)).^2 + (p(:,3)-recpar(it,3)).^2);
%     [~,a] = sort(rd);

    
%     index = a(1:N);
    
%     h     = rd(index(end));
%     Nind  = length(index);
%     wi    = diag(exp(-rd(index).^2/(h.^2)));
%     X     = ones(Nind,4);
%     X(:,1:3) = p(index,:);
    
    h     = max(dist(it,:));
    index = idx(it,:);
    wi    = diag(exp(-dist(it,:).^3/(h.^3)));
    X(:,1:3) = p(index,:);
    
    Rhs(:,:,it) = X'*wi;
    A1(:,:,it)  = X'*wi*X;
    
    jj(((it-1)*N+1):it*N)=index;
    
%     A=inv(X'*wi*X)*X'*wi;
    
%     jj(((it-1)*N+1):it*N)=index;
    
    
%     Lx(it,index) = A(1,:);
%     ABCD(it,:) =  (X'*wi*X)\(X'*wi*f(index));
    
end
% toc;

% tic;
A1inv = inv_4th_Order_Sym(A1);
% toc;

% tic;
for it = 1:Ndi
    A = A1inv(:,:,it)*Rhs(:,:,it);
    valx(((it-1)*N+1):it*N)=A(1,:);
    valy(((it-1)*N+1):it*N)=A(2,:);
    valz(((it-1)*N+1):it*N)=A(3,:);
end
% toc;

Lx = sparse(ii,jj,valx,Ndi,Nm);
Ly = sparse(ii,jj,valy,Ndi,Nm);
Lz = sparse(ii,jj,valz,Ndi,Nm);



