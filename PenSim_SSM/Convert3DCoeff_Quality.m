function [ Coeff_Quality_3D ] = Convert3DCoeff_Quality( Coeff_Quality, J_k, PastHorizon, ny, nu )
%convert the state space model and CVA to 3D ARX coeffcients
%   input:
%   Coeff_Quality: Q = x*Quality_Coeff
%   J_k: canonical loading
%   PastHorizon: a scalar indicating the lag
%   ny: dimension of y
%   nu: dimension of u
%   output:
%   Coeff_ARX_3D: MxPxPastHorizon
%   M=size(y)+size(u)

%   Q = size(q)

Coeff_temp = J_k'*Coeff_Quality;

[temp1, Q] = size(Coeff_temp);
M = temp1/PastHorizon;

if(M ~= ny+nu)
    error('ny + nu must equal to M');
end

Coeff_Quality_3D = zeros(M,Q,PastHorizon);

for i = 1:PastHorizon
    for j = 1:ny
        Coeff_Quality_3D(j,:,i) = Coeff_temp((PastHorizon-i)*ny+j,:);
        
    end
end

for i = 1:PastHorizon
    for j = ny+1:M
        Coeff_Quality_3D(j,:,i) = Coeff_temp(PastHorizon*ny+(PastHorizon-i)*nu+j-ny,:);
        
    end
end


end

