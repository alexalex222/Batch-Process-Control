function [ X ] = ScaleBack( Xscld, mean, sigma )
%Transform the scaled data to the raw data
%   input:
%   Xscld: the scaled data
%   mean: mean of the data
%   sigma: standard deviation
%   output
%   X: the raw data

X = Xscld.*repmat(sigma,size(Xscld,1),1)+repmat(mean, size(Xscld,1),1);


end

