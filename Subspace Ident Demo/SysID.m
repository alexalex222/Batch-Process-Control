function [ A_hat, B_hat, C_hat, D_hat, mt ] = SysID( P, F, Q, k, Yu, Uu  )
%Identify the linear system by using least squares
%   Input
%   P: past data (y,u)
%   F: future data (y)
%   k: dimension of etimated state
%   Yu: unfolded measurements
%   Uu: unfolded inputs
%   Output
%   A_hat: estimated A
%   B_hat: estimated B
%   C_hat: estimated C
%   D_hat: estimated D
%   Linear Time Invariant System
%   x(k+1)=Ax(k)+Bu(k)+w
%   y(k)=Cx(k)+Du(k)+Ew+v

[sample,~] = size(P);

[ ~, ~, ~, J_k ] = CVA( P, F, Q, k );


% Estimated States
% mt = (P-repmat(mean(P),sample,1))*J_k';
mt = P*J_k';

% Dimension of y
[~,ny]=size(Yu);
% Dimension of u
[~,nu]=size(Uu);
% Dimension of x
[~,nx]=size(mt);

XY = [mt(2:end,:),Yu(1:end-1,:)];
XU = [mt,Uu];

CovXYXU = (XY-repmat(mean(XY),sample-1,1))'*(XU(1:end-1,:)-repmat(mean(XU),sample-1,1))/(sample-2);
ABCD = CovXYXU*inv(cov(XU));

A_hat=ABCD(1:nx,1:nx);
B_hat=ABCD(1:nx,nx+1:nx+nu);
C_hat=ABCD(nx+1:nx+ny,1:nx);
D_hat=ABCD(nx+1:nx+ny,nx+1:nx+nu);

end

