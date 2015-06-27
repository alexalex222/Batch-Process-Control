%%%% gbngen.m; Matlab *.m file for generating a GBN signal
function U = gbngen(N,Tsw,Seed,umin,umax)
%
% This function generates GBN a signal U in [-1 1] with Tmin = 1 (samples
% to keep the signal constant)
%
% Input arguments:
% N : Number of Scimples
% Tsw : Average switching time in samples
% Seed : Seed for random generator
%
if nargin < 2
    error('Not enough input arguments, try again')
end
% Detemine switching probability
psw = 1/Tsw;

if nargin > 2
    rand('seed',Seed);
end;
R = rand(N,size(umin,1));

% Determine the initial value:
P_M(R(1,:)>.5) = 1;
P_M(R(1,:)<=.5) = -1;

U = zeros(N,size(umin,1));
for k=1:N
%     if (R(k) < psw)
%         P_M = -P_M;
%     end
    P_M(R(k,:) < psw) = -P_M(R(k,:)<psw);
    U(k,:) = P_M;
end

umin_all = repmat(umin',N,1);
umax_all = repmat(umax',N,1);
U(U == -1) = umin_all(U==-1);
U(U == 1) = umax_all(U==1);
U = U';

end