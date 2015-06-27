function[opts] = BatchSysOpts(varargin)
    % check for errors in input 
    if(mod(nargin,2) ~= 0)
        error('The number of inputs must be even');
    end

    % Initially set all values to their default.
    opts.nx = 4; % system size
    opts.nu = 3;
    opts.ny = 3;
    opts.nb = 100; % number of batches to genereate
    opts.IC_STD = 1.0; % initial condition standard deviation
    opts.meas_STD = 0; % measurement noise standard devaition
    opts.isAStable = true;
    opts.isRealEigenValues = true;
    opts.isZeroInitialState = false;
    opts.samt = 1;
    opts.batchLength = 100;
    opts.inputTrend = 'randomwalk'; % Can be 'stationary', 'randomwalk', 'linear'
    opts.inputVar = 'whitenoise'; % can be 'whitenoise', 'PRBS', or 'none'
    opts.whiteNoiseMag = 2; % magnitude for white noise in input variation
    opts.prbsMin = -1;
    opts.prbsMax = 1;
    opts.prbsTsw = 5;
    
    % Change values specified by user
    for i = 1:2:nargin
        switch varargin{i}
            case 'nx'
                CheckNumeric(varargin{i},varargin{i+1});
                opts.nx = varargin{i+1};
            case 'ny'
                CheckNumeric(varargin{i},varargin{i+1});
                opts.ny = varargin{i+1};
            case 'nu' 
                CheckNumeric(varargin{i},varargin{i+1});
                opts.nu = varargin{i+1};
            case 'nb' 
                CheckNumeric(varargin{i},varargin{i+1});
                opts.nb = varargin{i+1};
            case 'IC_STD'
                CheckNumeric(varargin{i},varargin{i+1});
                opts.IC_STD = varargin{i+1};
            case 'meas_STD'
                CheckNumeric(varargin{i},varargin{i+1});
                opts.meas_STD = varargin{i+1};
            case 'isAStable'
                opts.isAStable = varargin{i+1};
            case 'isRealEigenValues'
                opts.isRealEigenValues = varargin{i+1};
            case 'isZeroInitialState'
                opts.isZeroInitialState = varargin{i+1};
            case 'samt'
                CheckNumeric(varargin{i},varargin{i+1});
                opts.samt = varargin{i+1};
            case 'batchLength'
                CheckNumeric(varargin{i},varargin{i+1});
                opts.batchLength = varargin{i+1};
            case 'inputTrend'
                if(strcmp(varargin{i+1}, 'randomwalk') ||  ...
                    strcmp(varargin{i+1}, 'stationary') || ...
                    strcmp(varargin{i+1}, 'linear'))
                    opts.inputTrend = varargin{i+1};
                else
                    error('inputTrend must be ''randomwalk'', ''stationary'', or ''linear''')
                end
            case 'inputVar'
                if(strcmp(varargin{i+1}, 'whitenoise') ||  ...
                    strcmp(varargin{i+1}, 'PRBS') || ...
                    strcmp(varargin{i+1}, 'none'))
                    opts.inputVar = varargin{i+1};
                else
                    error('inputTrend must be ''whitenoise'', ''PRBS'', or ''none''')
                end
            case 'whiteNoiseMag'
                CheckNumeric(varargin{i},varargin{i+1});
                opts.whiteNoiseMag = varargin{i+1};
            case 'prbsMin'
                CheckNumeric(varargin{i},varargin{i+1});
                opts.prbsMin = varargin{i+1};
            case 'prbsMax'
                CheckNumeric(varargin{i},varargin{i+1});
                opts.prbsMax = varargin{i+1};
            case 'prbsTsw'
                CheckNumeric(varargin{i},varargin{i+1});
                opts.prbsTsw = varargin{i+1};                
            otherwise
                error([varargin{i} ' is not a recognized option']);
        end
    end
end

function[] = CheckNumeric(str, obj)
    if(numel(obj) ~= 1 || ~isnumeric(obj))
        error([str ' must be a scaler numeric value']);
    end
end
