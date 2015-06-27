function [Results] = Step_Forward(obj, varargin)
    %%%% Check status of model (ie check to see if current and model has
    %%%% been fit)  %%%%
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
    

    %%%% Unfold past data %%%%
    X = zeros(1,sum(sum(obj.Lags)));
    col = 1;
    for typidx = 1:size(obj.Types,2)
        for i = 1:size(obj.Measurements{typidx},2)
            varidx = obj.Measurements{typidx}(i);
            for lag = 1:obj.Lags(typidx,varidx)
                X(col) = past_data{typidx}(varidx, end-(lag-1));  % note (lag -1) becuase the data is already past data
                col = col + 1;                                        %   ie: lag of 1 is the last value in the data set
            end
        end
    end
    
    %%%% Scale and center past data %%%%
    X_scld = (X- obj.Measure_Means)./obj.Measure_Std;
    
    %%%% Determine Membership %%%%
    D = Mult_Mod.distfcm(obj.Centers, X_scld);
    tmp = D.^(-2);
    MF = tmp./(ones(size(obj.Centers,1), 1)*sum(tmp));
    
    %%%% Apply the model %%%%
    tot_measured = sum(sum(obj.Lags));
    phi = reshape(repmat(MF,1,(tot_measured+1))',1,(tot_measured+1)*obj.Opt_Clusters)...
        .*repmat([1 X_scld],1,obj.Opt_Clusters);
    
    if(strcmpi(obj.Regression_Type, 'PLS'))
%        Y_scld = [1 phi]*obj.Coeffs;
         Y_scld = [phi]*obj.Coeffs;
    else
        Y_scld = phi*obj.Coeffs;
    end
    
    %%%% Un-scale result and fill output structure %%%%
    Y = Y_scld.*obj.Predicted_Std + obj.Predicted_Means;
    
    Results = cell(size(obj.Types,2),1);
    
    k = 1;
    for typidx = 1:size(obj.Types,2)
       for i = 1:size(obj.Predicted{typidx},2)
          varidx = obj.Predicted{typidx}(i);
          Results{typidx}(varidx,1) = Y(k);
          k = k+1;
       end
    end
    
end