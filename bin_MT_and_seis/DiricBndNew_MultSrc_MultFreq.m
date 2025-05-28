function [K, b]   = DiricBndNew_MultSrc_MultFreq(K,b,Bdirc,freq,nfreq,ns)

%% Add Dirichlet boundary condition which is stored in BGdirc in to Kglob and bglob

if isempty(Bdirc)
    fprintf('No boundary condition is added\n');
    return;
end

for ifreq = 1:nfreq
    for isrc = 1:ns
        b{ifreq,isrc}             = b{ifreq,isrc} - (real(K(:,Bdirc(:,1)))+1i*(freq(ifreq)/(freq(1)))*(K(:,Bdirc(:,1))))*Bdirc(:,2);
        b{ifreq,isrc}(Bdirc(:,1)) = Bdirc(:,2);
    end
end




K(Bdirc(:,1),:)   = 0;
K(:,Bdirc(:,1))   = 0;

nk = size(K,1);
ii = Bdirc(:,1);
Kbnd = sparse(ii,ii,ones(size(ii)),nk,nk);

K = K+Kbnd;