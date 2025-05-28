function [x,flag,relres,iter,resvec,lsvec] = lsqr_Implicit_A_MaxIter(A1_dense,A2_sparse,b,tol,maxit)
%Modified by Hongzhu Cai
% A modified version of MATLAB lsqr
%assume input is matrix
%input matrix A has two part A=[A1_dense,A2_sparse], the two section is
%seperated for memroy issue, dense part is A1_dense, the sparse part is
%A2_sparse.


n = size(A1_dense,2);
m = size(A1_dense,1) + size(A2_sparse,1);

nKnown = true;
xInit = false;    
x0 = zeros(size(b));




% Set up for the method
n2b = norm(b);                     % Norm of rhs vector, b
flag = 1;
tolb = tol * n2b;                  % Relative tolerance
u = b;
beta = norm(u);
% Norm of residual r=b-A*x is estimated well by prod_i abs(sin_i)
normr = beta;
if beta ~= 0
    u = u / beta;
end
c = 1;
s = 0;
phibar = beta;
% v = A'*u;
v = Ax_Implicit(A1_dense,A2_sparse,u,1);
if ~nKnown
    % How many columns does A have?  Same as entries of v = A'*u.
    n = size(v,1);
end
% And make sure we have x if it was not initialized.
if ~xInit
    x = zeros(n,1);
end

alpha = norm(v);
if alpha ~= 0
    v = v / alpha;
end
d = zeros(n,1);

% norm((A*inv(M))'*r) = alpha_i * abs(sin_i * phi_i)
normar = alpha * beta;

% Check for exact solution
if (normar == 0)               % if alpha_1 == 0 | beta_1 == 0
    flag = 0;                  % a valid solution has been obtained
    relres = normr/n2b;
    iter = 0;                  % no iterations need be performed
    resvec = normr;            % resvec(1) = norm(b-A*x) = norm(0)
    lsvec = zeros(0,1);        % no estimate for norm(A*inv(M),'fro') yet
    return
end



% Poorly estimate norm(A*inv(M),'fro') by norm(B_{ii+1,ii},'fro')
% which is in turn estimated very well by
% sqrt(sum_i (alpha_i^2 + beta_{ii+1}^2))
norma = 0;
% norm(inv(A*inv(M)),'fro') = norm(D,'fro')
% which is poorly estimated by sqrt(sum_i norm(d_i)^2)
sumnormd2 = 0;
resvec = zeros(maxit+1,1);     % Preallocate vector for norm of residuals
resvec(1) = normr;             % resvec(1,1) = norm(b-A*x0)
lsvec = zeros(maxit,1);        % Preallocate vector for least squares estimates
stag = 0;                      % stagnation of the method
iter = maxit;                  % Assume lack of convergence until it happens
maxstagsteps = 3;

% loop over maxit iterations (unless convergence or failure)

for ii = 1 : maxit
   z = v;
   
    Az = Ax_Implicit(A1_dense,A2_sparse,z,0);
    u = Az - alpha * u;
%     u = A*z - alpha * u;
    beta = norm(u);
    u = u / beta;
    norma = norm([norma alpha beta]);
    lsvec(ii) = normar / norma;
    thet = - s * alpha;
    rhot = c * alpha;
    rho = sqrt(rhot^2 + beta^2);
    c = rhot / rho;
    s = - beta / rho;
    phi = c * phibar;
    if (phi == 0)              % stagnation of the method
        stag = 1;
    end
    phibar = s * phibar;
    d = (z - thet * d) / rho;
    sumnormd2 = sumnormd2 + (norm(d))^2;
    
    % Check for stagnation of the method
    if abs(phi)*norm(d) < eps*norm(x)
        stag = stag + 1;
    else
        stag = 0;
    end
    
    if normar/(norma*normr) <= tol % check for convergence in min{|b-A*x|}
        flag = 0;
        iter = ii-1;
        resvec = resvec(1:iter+1);
        lsvec = lsvec(1:iter+1);
%         break
    end
    
    if normr <= tolb           % check for convergence in A*x=b
        flag = 0;
        iter = ii-1;
        resvec = resvec(1:iter+1);
        lsvec = lsvec(1:iter+1);
        break
    end
    
    if stag >= maxstagsteps
        flag = 3;
        iter = ii-1;
        resvec = resvec(1:iter+1);
        lsvec = lsvec(1:iter+1);
%         break
    end
    
    x = x + phi * d;
    normr = abs(s) * normr;
    resvec(ii+1) = normr;
    vt = Ax_Implicit(A1_dense,A2_sparse,u,1);
%     vt = A'*u;
    v = vt - beta * v;
    alpha = norm(v);
    v = v / alpha;
    normar = alpha * abs( s * phi);
    
end                            % for ii = 1 : maxit

if flag == 1
    if normar/(norma*normr) <= tol % check for convergence in min{|b-A*x|}
        flag = 0;
        iter = maxit;
    end
    
    if normr <= tolb           % check for convergence in A*x=b
        flag = 0;
        iter = maxit;
    end
end

relres = normr/n2b;
disp(['Lsqr runs ' num2str(ii) 'iterations with residual of ' num2str(relres)]);

function Ax = Ax_Implicit(A1_dense,A2_sparse,x,trans)

n1 = size(A1_dense,1);
% n2 = size(A2_sparse,1);

if trans
    Ax = A1_dense'*x(1:n1) + A2_sparse'*x((n1+1):end);
else
    A1x = A1_dense*x ;
    A2x = A2_sparse*x;
    Ax  = [A1x;A2x];
end