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
R = rand(N,1);

% Determine the initial value:
if R(1) > 0.5,
    P_M = 1;
else
    P_M = -1;
end;

U = zeros(N,1);
for k=1:N
    if (R(k) < psw)
        P_M = -P_M;
    end
    U(k) = P_M;
end

U(U == -1) = umin;
U(U == 1) = umax;
U = U';

end