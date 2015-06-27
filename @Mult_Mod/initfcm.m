function U = initfcm(cluster_n, data_n)
%INITFCM Generate initial fuzzy partition matrix for fuzzy c-means clustering.
%	U = INITFCM(CLUSTER_N, DATA_N) randomly generates a fuzzy partition
%	matrix U that is CLUSTER_N by DATA_N, where CLUSTER_N is number of
%	clusters and DATA_N is number of data points. The summation of each
%	column of the generated U is equal to unity, as required by fuzzy
%	c-means clustering.
%
%       See also DISTFCM, FCMDEMO, IRISFCM, STEPFCM, and FCM.

%	Roger Jang, 12-1-94.
%       Copyright (c) 1994-95  by The MathWorks, Inc.
%       $Revision: 1.4 $  $Date: 1995/02/17 13:08:10 $

U = rand(cluster_n, data_n);
col_sum = sum(U);
U = U./col_sum(ones(cluster_n, 1), :);
