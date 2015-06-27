clear variables;
close all;

%__________________________________________________________________________
% User parameters
SeedBase = 100; % use this to try different systems
nx = 9; % system size
nu = 3;
ny = 3;
nb = 100; % number of batches to genereate
InitialConditionSTD = 1.0;
MeasurementNoiseSTD = 0;
isAStable = true;
isRealEigenValues = true;
isZeroStateInitial = false;
samt = 1;
batchLength = 100;
inputTrend = 'linear'; % Can be 'stationary', 'randomwalk', 'linear'
inputVariation = 'whitenoise'; % can be 'whitenoise', 'PRBS', or 'none'
whiteNoiseMagnitude = 2;
prbsMin = -1;
prbsMax = 1;
prbsTsw = 5;

%==========================================================================
% * * * DO NOT MODIFY CODE BELOW THIS LINE WITHOUT DISCUSSION * * *
%__________________________________________________________________________


% Generate system
% Create random eigen vectors
rand('seed', SeedBase);
v = orth(-1+2*rand(nx));

% select eigen values
rand('seed',2*SeedBase);
if(isRealEigenValues)
    if(isAStable)
        angles = repmat(pi,nx,1);
    else
        angles = 0 + pi* round(rand(nx,1));
    end
else
    if(isAStable)
        angles = pi/2 + pi*rand(nx,1);
    else
        angles = -pi + 2*pi*rand(nx,1);
    end
end
rads = 10*rand(nx,1);

if(isRealEigenValues)
    lambda = cos(angles).*rads;
else
    lambda = cos(angles).*rads + sin(angles).*rads*i;
end

% generate A matrix
Acont = v*diag(lambda)*v';

% select a B matrix
rand('seed', 3*SeedBase)
Bcont = -.5 + 1*rand(nx,nu);

% select a C matrix
Ccont = -2 + 4*rand(ny, nx);

% Check controlability and observability


% Get discrete matricies 
sys = ss(Acont, Bcont, Ccont, zeros(ny, nu));

sys = c2d(sys, samt);

A = sys.a;
B = sys.b;
C = sys.c;

% Generate data
bdb = BatchDatabase(nx,ny,nu, 0, batchLength, 1);


bdb.AddICStd(repmat(InitialConditionSTD, nx,1));

% Set up measurement noise
bdb.AddMeasStd(repmat(MeasurementNoiseSTD, ny, 1));

% generate nominal batch
x = zeros(nx,batchLength+1);
y = zeros(ny, batchLength+1);
u = zeros(nu, batchLength);
q = zeros(0,batchLength+1);

% generate initial state
if(isZeroStateInitial)
   x(:,1) = zeros(nx,1);
else
   rand('seed', 4*SeedBase);
   x(:,1) = -3 + 6*rand(nx,1);
end
y(:,1) = C*x(:,1);

% generate nominal input (this is just the tend)
switch(inputTrend)
    case 'stationary'
        u = zeros(nu,batchLength);
        
    case 'randomwalk'
        rand('seed', 5*SeedBase);
        u(:,1) = -3 + 6*rand(nu,1);
        for i = 1:batchLength-1
           u(:,i+1) = u(:,i) +(-1+2*rand(nu,1));
        end
        
    case 'linear'
        rand('seed', 6*SeedBase);
        m = -.1 + .2*rand(nu,1);
        u(:,1) = -3 + 6*rand(nu,1);
        u(:,2:end) = repmat(u(:,1),1,batchLength-1) + repmat(m,1,batchLength-1).*repmat(2:batchLength,nu,1);
end

% generate nominal data
for i = 1:batchLength
   x(:,i+1) = A*x(:,i) + B*u(:,i);
   y(:,i+1) = C*x(:,1+i);
end
bdb.AddNomBatch(x,y,u,q);

% get all initial conditions
X0 = bdb.GetICs(nb,7*SeedBase);

% generate data
for b = 1:nb
    x = zeros(nx,batchLength+1);
    y = zeros(ny, batchLength+1);
    u = bdb.data.unom;
    q = zeros(0,batchLength+1);
    
    switch(inputVariation)
        case 'none'
            % do nothing
        case 'whitenoise'
            u = u+ (-whiteNoiseMagnitude +2*whiteNoiseMagnitude*rand(nu, batchLength));
        case 'PRBS'
            u = u + gbngen(batchLength, prbsTsw, 8*SeedBase*b, repmat(prbsMin,nu,1), repmat(prbsMax,nu,1));
    end
    
    % get measurement noise 
    noise = bdb.GenNoise(b*123);
    
    x(:,1) = X0(b,:)';
    y(:,1) = C*x(:,1) + noise(:,1);
    
    for i = 1:batchLength
        x(:,i+1) = A*x(:,i) + B*u(:,i);
        y(:,i+1) = C*x(:,1+i) + noise(:,i+1);
    end
    
    bdb.AddBatch(x,y,u,q);
end

% Plot results
plotmtxdim = ceil(sqrt(nx));

for i = 1:nx
    subplot(plotmtxdim, plotmtxdim, i);
    bdb.PlotTraj('state', i, 1:bdb.nb,'r');
    hold on;
    plot(0:batchLength,bdb.data.xnom(i, :),'k');
    title(['State ' num2str(i)])
end

figure();
plotmtxdim = ceil(sqrt(ny));

for i = 1:ny
    subplot(plotmtxdim, plotmtxdim, i);
    bdb.PlotTraj('output', i, 1:bdb.nb,'r');
    hold on;
    plot(0:batchLength,bdb.data.ynom(i, :),'k');
    title(['Output ' num2str(i)])
end

figure();
plotmtxdim = ceil(sqrt(nu));

for i = 1:ny
    subplot(plotmtxdim, plotmtxdim, i);
    bdb.PlotTraj('input', i, 1:bdb.nb,'r');
    hold on;
    stairs(0:batchLength-1,bdb.data.unom(i, :),'k');
    title(['Input ' num2str(i)])
end

