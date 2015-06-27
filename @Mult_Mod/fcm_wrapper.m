function [MF C centers J]  = fcm_wrapper(X,L,nFCMiter)

N = size(X,1);

% populate matrices:
centers = zeros(L,size(X,2),nFCMiter); % FCM centers
U = zeros(L,N,nFCMiter);               % FCM memberships
J = zeros(nFCMiter,1);                 % FCM objective function

for i=1:nFCMiter
    [c_curr, U_curr, J_curr] = Mult_Mod.fcm(X,L,[NaN 4000000 1e-8 0]); % call FCM with options (NaN selects default)
    centers(:,:,i) = c_curr;
    U(:,:,i) = U_curr;
    J(i,1) = J_curr(end);
    fprintf('--- FCM iteration number: %d of %d \n', i, nFCMiter);
end
[dummy Jmin_idx] = min(J); %#ok<ASGLU>
MF = U(:,:,Jmin_idx);
C = centers(:,:,Jmin_idx);
end