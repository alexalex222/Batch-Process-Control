function [ Coeff_ARX_3D ] = Convert3DCoeff_ARX( C, J_k, PastHorizon )
%convert the state space model and CVA to 3D ARX coeffcients
%   input:
%   C: y = Cx
%   J_k: canonical loading
%   PastHorizon: a scalar indicating the lag
%   output:
%   Coeff_ARX_3D: MxPxPastHorizon
%   M=size(y)+size(u)
%   P=size(y)

Coeff_temp = (C*J_k)';

[temp1, P] = size(Coeff_temp);
M = temp1/PastHorizon;
Nu = temp1/PastHorizon - P;

Coeff_ARX_3D = zeros(M,P,PastHorizon);

for i = 1:PastHorizon
    for j = 1:P
        Coeff_ARX_3D(j,:,i) = Coeff_temp((PastHorizon-i)*P+j,:);
        
    end
end

for i = 1:PastHorizon
    for j = P+1:M
        Coeff_ARX_3D(j,:,i) = Coeff_temp(PastHorizon*P+(PastHorizon-i)*Nu+j-P,:);
        
    end
end


end

