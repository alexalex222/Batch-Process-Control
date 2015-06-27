function [SPE Tsqrd] = getSPEandTsqrd(obj, varargin)

if(~obj.Model_Fitted)
       error('No model has been fitted, use the ''Fit_Model'' method after model structure is completed.');
    end
    if(~obj.Model_Current)
        response = input(['Changes may have been'...
            ' made to parameters of the model since the last time the model was fit.\n'...
            ' Do you still want to use the model? Y/N [N] '],'s');
        if(isempty(response))
            response = 'N';
        end
        if(strcmpi(response, 'N'))
            error('Model not current');
        end
    end

    %%%% Get data that was passed in %%%%
    if(mod(size(varargin,2),2) ~=0)
        error('There is a problem with the passed variables');
    end
    
    if(strcmp(varargin{1},'struct'))
        past_data = varargin{2};
    else
        nvarpass = size(varargin,2)/2;
        past_data = cell(size(obj.Types,2),1);
        for n = 1:nvarpass
            type = varargin{(n-1)*2+1};
            typidx = getTypeIndex(obj, type);
            past_data{typidx} = varargin{2*n};
        end

        % Check that all needed data types were passed
        for typidx = 1:size(obj.Types,2)
            if(isempty(past_data{typidx}) && ~isempty(obj.Measurements{typidx}));
                error(['The values for past ''' obj.Types{typidx} ''' were not specified']);
            end
        end
    end
    
    %%%% determine the number of times the model could have been applied
    numTimes = intmax;
    for i = 1:size(obj.Types, 2)
        if(~isempty(obj.Measurements{i}))
            if(size(past_data{i},2)-(max(obj.Lags(i,:))-1)<numTimes)
                numTimes = size(past_data{i},2)-(max(obj.Lags(i,:))-1);
            end
        end
    end

    %%%% Unfold past data %%%%
    X = zeros(numTimes,sum(sum(obj.Lags)));
    
    for t = 1:numTimes
        col = 1;
        for typidx = 1:size(obj.Types,2)
            for i = 1:size(obj.Measurements{typidx},2)
                varidx = obj.Measurements{typidx}(i);
                for lag = 1:obj.Lags(typidx,varidx)
                    X(t,col) = past_data{typidx}(varidx, t+(lag-1));  % note (lag -1) becuase the data is already past data
                    col = col + 1;                                        %   ie: lag of 1 is the last value in the data set
                end
            end
        end
    end
    
    %%%% Scale and center past data %%%%
    X_scld = (X- repmat(obj.Measure_Means,numTimes,1))./repmat(obj.Measure_Std,numTimes,1);
    
    for t = 1:numTimes
        %%%% Determine Membership %%%%
        D = Mult_Mod.distfcm(obj.Centers, X_scld(t,:));
        tmp = D.^(-2);
        MF = tmp./(ones(size(obj.Centers,1), 1)*sum(tmp));

        %%%% Apply the model %%%%
        tot_measured = sum(sum(obj.Lags));
        phi = reshape(repmat(MF,1,(tot_measured+1))',1,(tot_measured+1)*obj.Opt_Clusters)...
            .*repmat([1 X_scld(t,:)],1,obj.Opt_Clusters);

        xs = phi*obj.wstar;

        xhat = xs*obj.loading';
        err = phi- xhat;
        SPE(t) = sum(err.^2,2);

        Tsqrd(t) = sum((xs./obj.sa).^2);
    end