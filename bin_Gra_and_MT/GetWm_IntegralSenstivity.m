function Wm = GetWm_IntegralSenstivity_p(J,p)
Nm = size(J,2);
Wm_diag = zeros(Nm,1);
for nm=1:Nm
    Jn = J(:,nm)*p(Nm);
    Wm_diag(nm)=(Jn'*Jn)^1.0;
  % Wm_diag(nm)=(Jn'*Jn)^0.5;
%Wm_diag(nm)=((Jn'*Jn))^(0.22);
end
% 
 Wm = spdiags(Wm_diag,0,Nm,Nm);

% function Wm = GetWm_IntegralSenstivity(J)
% A = diag(J'*J);
% B = A.^(0.25);
% Wm = diag(B);

% Nm = size(J,2);
% Wm_diag = zeros(Nm,1);
% for nm=1:Nm
%     Wm_diag(nm)=((J(:,nm))'*J(:,nm))^0.25;
% end
% 
% Wm = spdiags(Wm_diag,0,Nm,Nm);