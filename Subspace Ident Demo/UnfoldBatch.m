function [ P, F, Q, Xu, Yu, Uu,nb ] = UnfoldBatch( X, Y, U, PastHorizon, FutureHorizon )
%Unfold and scale the batch data for subspace identification
%   input
%   X: multiway data for states
%   Y: multiway data for measurements 
%   U: multiway data for control inputs
%   PastHorizon:  the lag of past data
%   FutureHorizon: the lag of fuuture data
%   Output
%   P: stacked past input/ouputs
%   F: stacked future outputs
%   Q: stacked future inputs
%   Xu: unfolded state
%   Yu: unfolded measurement
%   Uu: unfolded control input

% delete the inital conditions in X and Y
Y = Y(:,2:end,:);
U = U;
X = X(:,2:end,:);

% The number of batches in the multiway data
[nx,~,nb]=size(X);
ny = size(Y,1);
[nu,nFE,~] = size(U);

% For each batch, the available data points from PastHorizon+1 to
% N-FutureHorizon+1

% P = zeros(nb*(bdb.nFE-PastHorizon-FutureHorizon+1),(bdb.ny+bdb.nu)*PastHorizon);
% F = zeros(nb*(bdb.nFE-PastHorizon-FutureHorizon+1),bdb.ny*FutureHorizon);
% Xu = zeros(nb*(bdb.nFE-PastHorizon-FutureHorizon+1),bdb.nx);
% Yu = zeros(nb*(bdb.nFE-PastHorizon-FutureHorizon+1),bdb.ny);
% Uu = zeros(nb*(bdb.nFE-PastHorizon-FutureHorizon+1),bdb.nu);
P = zeros(nb*(nFE-PastHorizon-FutureHorizon),(ny+nu)*PastHorizon);
F = zeros(nb*(nFE-PastHorizon-FutureHorizon),ny*FutureHorizon);
Q = zeros(nb*(nFE-PastHorizon-FutureHorizon),nu*FutureHorizon);
Xu = zeros(nb*(nFE-PastHorizon-FutureHorizon),nx);
Yu = zeros(nb*(nFE-PastHorizon-FutureHorizon),ny);
Uu = zeros(nb*(nFE-PastHorizon-FutureHorizon),nu);
p = zeros(1,(ny+nu)*PastHorizon);
f = zeros(1,ny*FutureHorizon);
q = zeros(1,nu*FutureHorizon);

% In generating the batch data, u(1) is put into ODEs to get x(1).
% Therefore, u(1) is u0

for j=1:nb
    for i=1:nFE-PastHorizon-FutureHorizon
        tempyp = Y(:,i:i+PastHorizon-1,j)';
        tempyp = reshape(tempyp',1,size(tempyp,1)*size(tempyp,2));
        % Note that there is a shift 
        tempup = U(:,i+1:i+PastHorizon,j)';
        tempup = reshape(tempup',1,size(tempup,1)*size(tempup,2));
        p=[tempyp,tempup];
        P((j-1)*(nFE-PastHorizon-FutureHorizon+1)+i,:)=p;
        tempyf = Y(:,i+PastHorizon: i+PastHorizon+FutureHorizon-1,j)';
        tempyf = reshape(tempyf',1,size(tempyf,1)*size(tempyf,2));
        tempuf = U(:,i+PastHorizon+1: i+PastHorizon+FutureHorizon,j)';
        tempuf = reshape(tempuf',1,size(tempuf,1)*size(tempuf,2));
        f = tempyf;
        F((j-1)*(nFE-PastHorizon-FutureHorizon+1)+i,:)=f;
        q = tempuf;
        Q((j-1)*(nFE-PastHorizon-FutureHorizon+1)+i,:)=q;
        
        tempx = X(:,i+PastHorizon,j);
        Xu((j-1)*(nFE-PastHorizon-FutureHorizon+1)+i,:)=tempx;
        tempu = U(:,i+PastHorizon+1,j);
        Uu((j-1)*(nFE-PastHorizon-FutureHorizon+1)+i,:)=tempu;
        tempy = Y(:,i+PastHorizon,j);
        Yu((j-1)*(nFE-PastHorizon-FutureHorizon+1)+i,:)=tempy;
        
    end
end


end

