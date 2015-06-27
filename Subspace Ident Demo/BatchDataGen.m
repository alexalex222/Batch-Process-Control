function[bdb, LinearSystem] = BatchDataGen(opts, seed, plot_result)
%     %__________________________________________________________________________
%     % User parameters
%     SeedBase = 100; % use this to try different systems
%     opts.nx = 9; % system size
%     opts.nu = 3;
%     opts.ny = 3;
%     opts.nb = 100; % opts.number of batches to genereate
%     opts.IC_STD = 1.0;
%     opts.meas_STD = 0;
%     opts.isAStable = true;
%     opts.isRealEigenValues = true;
%     opts.isZeroInitialState = false;
%     opts.samt = 1;
%     opts.batchLength = 100;
%     opts.inputTrend = 'linear'; % Can be 'stationary', 'randomwalk', 'linear'
%     opts.inputVar = 'whitenoise'; % can be 'whitenoise', 'PRBS', or 'none'
%     opts.whiteNoiseMag = 2;
%     opts.prbsMin = -1;
%     opts.prbsMax = 1;
%     opts.prbsTsw = 5;

    %==========================================================================
    % * * * DO NOT MODIFY CODE BELOW THIS LINE WITHOUT DISCUSSION * * *
    %__________________________________________________________________________

    SeedBase = seed;

    % Generate system
    % Create random eigen vectors
    rand('seed', SeedBase);
    v = orth(-1+2*rand(opts.nx));

    % select eigen values
    rand('seed',2*SeedBase);
    if(opts.isRealEigenValues)
        if(opts.isAStable)
            angles = repmat(pi,opts.nx,1);
        else
            angles = 0 + pi* round(rand(opts.nx,1));
        end
    else
        if(opts.isAStable)
            angles = pi/2 + pi*rand(opts.nx,1);
        else
            angles = -pi + 2*pi*rand(opts.nx,1);
        end
    end
    rads = 10*rand(opts.nx,1);

    if(opts.isRealEigenValues)
        lambda = cos(angles).*rads;
    else
        lambda = cos(angles).*rads + sin(angles).*rads*j;
    end

    % generate A matrix
    Acont = v*diag(lambda)*v';

    % select a B matrix
    rand('seed', 3*SeedBase)
    Bcont = -.5 + 1*rand(opts.nx,opts.nu);

    % select a C matrix
    Ccont = -2 + 4*rand(opts.ny, opts.nx);

    % Check controlability and observability
    while rank(obsv(Acont,Ccont)) < opts.nx || rank(ctrb(Acont,Bcont)) < opts.nx
        SeedBase = SeedBase+1;
        rand('seed', 3*SeedBase);
        Bcont = -.5 + 1*rand(opts.nx,opts.nu);
        Ccont = -2 + 4*rand(opts.ny, opts.nx);
    end


    % Get discrete matricies 
    sys = ss(Acont, Bcont, Ccont, zeros(opts.ny, opts.nu));

    sys = c2d(sys, opts.samt);

    A = sys.a;
    B = sys.b;
    C = sys.c;
    
    LinearSystem.Acont = Acont;
    LinearSystem.Bcont = Bcont;
    LinearSystem.Ccont = Ccont;
    LinearSystem.A = A;
    LinearSystem.B = B;
    LinearSystem.C = C;
    
    
   

    % Generate data
    bdb = BatchDatabase(opts.nx,opts.ny,opts.nu, 0, opts.batchLength, 1);
    
    


    bdb.AddICStd(repmat(opts.IC_STD, opts.nx,1));

    % Set up measurement noise
    bdb.AddMeasStd(repmat(opts.meas_STD, opts.ny, 1));

    % generate nominal batch
    x = zeros(opts.nx,opts.batchLength+1);
    y = zeros(opts.ny, opts.batchLength+1);
    u = zeros(opts.nu, opts.batchLength);
    q = zeros(0,opts.batchLength+1);

    % generate initial state
    if(opts.isZeroInitialState)
       x(:,1) = zeros(opts.nx,1);
    else
       rand('seed', 4*SeedBase);
       x(:,1) = -3 + 6*rand(opts.nx,1);
    end
    y(:,1) = C*x(:,1);

    % generate nominal input (this is just the tend)
    switch(opts.inputTrend)
        case 'stationary'
            u = zeros(opts.nu,opts.batchLength);

        case 'randomwalk'
            rand('seed', 5*SeedBase);
            u(:,1) = -3 + 6*rand(opts.nu,1);
            for i = 1:opts.batchLength-1
               u(:,i+1) = u(:,i) +(-1+2*rand(opts.nu,1));
            end

        case 'linear'
            rand('seed', 6*SeedBase);
            m = -.1 + .2*rand(opts.nu,1);
            u(:,1) = -3 + 6*rand(opts.nu,1);
            u(:,2:end) = repmat(u(:,1),1,opts.batchLength-1) + repmat(m,1,opts.batchLength-1).*repmat(2:opts.batchLength,opts.nu,1);
    end

    % generate nominal data
    for i = 1:opts.batchLength
       x(:,i+1) = A*x(:,i) + B*u(:,i);
       y(:,i+1) = C*x(:,1+i);
    end
    bdb.AddNomBatch(x,y,u,q);

    % get all initial conditions
    X0 = bdb.GetICs(opts.nb,7*SeedBase);

    % generate data
    for b = 1:opts.nb
        x = zeros(opts.nx,opts.batchLength+1);
        y = zeros(opts.ny, opts.batchLength+1);
        u = bdb.data.unom;
        q = zeros(0,opts.batchLength+1);

        switch(opts.inputVar)
            case 'none'
                % do nothing
            case 'whitenoise'
                u = u+ (-opts.whiteNoiseMag +2*opts.whiteNoiseMag*rand(opts.nu, opts.batchLength));
            case 'PRBS'
                u = u + gbngen(opts.batchLength, opts.prbsTsw, 8*SeedBase*b, repmat(opts.prbsMin,opts.nu,1), repmat(opts.prbsMax,opts.nu,1));
        end

        % get measurement noise 
        noise = bdb.GenNoise(b*123);

        x(:,1) = X0(b,:)';
        y(:,1) = C*x(:,1) + noise(:,1);

        for i = 1:opts.batchLength
            x(:,i+1) = A*x(:,i) + B*u(:,i);
            y(:,i+1) = C*x(:,1+i) + noise(:,i+1);
        end

        bdb.AddBatch(x,y,u,q);
    end

    if(plot_result)
        % Plot results
        plotmtxdim = ceil(sqrt(opts.nx));

        for i = 1:opts.nx
            subplot(plotmtxdim, plotmtxdim, i);
            bdb.PlotTraj('state', i, 1:bdb.nb,'r');
            hold on;
            plot(0:opts.batchLength,bdb.data.xnom(i, :),'k');
            title(['State ' num2str(i)])
        end

        figure();
        plotmtxdim = ceil(sqrt(opts.ny));

        for i = 1:opts.ny
            subplot(plotmtxdim, plotmtxdim, i);
            bdb.PlotTraj('output', i, 1:bdb.nb,'r');
            hold on;
            plot(0:opts.batchLength,bdb.data.ynom(i, :),'k');
            title(['Output ' num2str(i)])
        end

        figure();
        plotmtxdim = ceil(sqrt(opts.nu));

        for i = 1:opts.nu
            subplot(plotmtxdim, plotmtxdim, i);
            bdb.PlotTraj('input', i, 1:bdb.nb,'r');
            hold on;
            stairs(0:opts.batchLength-1,bdb.data.unom(i, :),'k');
            title(['Input ' num2str(i)])
        end
    end
end