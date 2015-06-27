function [ J, L, D, J_k, SigmaXX, SigmaYY, Omega ] = CVA( P, F, Q, k )
%Canonical Variate Analysis
%   Input
%   X: past input/outputs
%   Y: future outputs
%   Q; future inputs
%   k: the dimension of subspace
%   Output
%   J: loading matrix for X
%   L: loading matrix for Y
%   D: diagonal matrix containing canonical correaltion

[sample,dx] = size(P);
% [~,dy] = size(Y);


% Future data correction
Beta_ARX = ols([P,Q], F);

% [~,~,~,~,Beta_ARX] = plsregress([P,Q], F, size([P,Q],2));
Omega = Beta_ARX(end-size(Q,2)+1:end,:);
% F = F - Q*Omega;

% SigmaXX = cov(X);
SigmaXX = (P-repmat(mean(P),sample,1))'*(P-repmat(mean(P),sample,1))/(sample-1);
SigmaYY = cov(F);
SigmaXY = (P-repmat(mean(P),sample,1))'*(F-repmat(mean(F),sample,1))/(sample-1);

SigmaXX = nearestSPD(SigmaXX);
SigmaYY = nearestSPD(SigmaYY);

[U, D, V] = svd(SigmaXX^(-0.5)*SigmaXY*SigmaYY^(-0.5));

J = U'*SigmaXX^(-0.5);
J_k = U(:,1:k)'*SigmaXX^(-0.5);
L = V'*SigmaYY^(-0.5);


end

