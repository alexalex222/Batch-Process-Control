function [ Beta ] = ols( X, Y )
%Ordinary least squares
%   Note that X and Y are raw data
%   Input
%   X: regressor data
%   Y: response data
%   Output:
%   Beta: least-squares estimation of regression parameter

Beta = inv(X'*X)*X'*Y;


end

