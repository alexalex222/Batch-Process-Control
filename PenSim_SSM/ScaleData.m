function [ Xscld ] = ScaleData( X, mean, sigma )
%Scale the data: Xscld = (X - mean)/sigma
%   input:
%   X: raw data (sample x variable)
%   mean: the mean of data
%   sigma: the standard deviation
%   output:
%   Xscld: the scaled data

Xscld = (X - repmat(mean,size(X,1),1))./repmat(sigma,size(X,1),1);


end

