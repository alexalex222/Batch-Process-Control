function [predicted_values] = Predict_Forward(obj, past_length, steps, varargin)
    
    %%%% Check model status %%%%
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
    
    % NOTE: passed data should contain the past data needed to initialize
    % the model and possibly some future data (ie control inputs)
    if(mod(size(varargin,2),2) ~=0)
        error('There is a problem with the passed variables');
    end
    
    if(strcmp(varargin{1},'struct'))
        user_data = varargin{2};
    else
        nvarpass = size(varargin,2)/2;
        user_data = cell(size(obj.Types,2),1);
        for n = 1:nvarpass
            type = varargin{(n-1)*2+1};
            typidx = getTypeIndex(obj, type);
            user_data{typidx} = varargin{2*n};
        end
    end
    
    for typidx = 1:size(obj.Types,2)
        if(isempty(user_data{typidx}) && ~isempty(obj.Measurements{typidx}));
            error(['The values for past ''' obj.Types{typidx} ''' were not specified']);
        end
    end
    
    %%%% Sort passed information into past and future data %%%%
    past_data = cell(size(obj.Types,2),1);
    future_data = cell(size(obj.Types,2),1);
    
    for typidx = 1:size(obj.Types,2)
        % for all variables used as measurements in the model
        for i = 1:size(obj.Measurements{typidx},2)
            varidx = obj.Measurements{typidx}(i);
            
            % Copy all past values from user data 
            past_data{typidx}(varidx,:) = user_data{typidx}(varidx,1:past_length);
            
            % For those values that are not predicted by the model, copy
            % future values from user data.  NOTE: In the case where an
            % output is not predicted by the model this is equivilant to saying that an
            % external perfect model is used in this prediction.  
            if(~ismember(varidx, obj.Predicted{typidx}))
                future_data{typidx}(varidx,:) = user_data{typidx}(varidx,past_length+1:end);
            end
        end
    end
    
    
    %%%% Predict future values %%%%
    predicted_values = cell(size(obj.Types,2),1);
    
    
    for k = past_length + 1:past_length+ steps
        % clear measured_values at each time step
        measured_values = cell(size(obj.Types,2),1);
        
        % compile all measurements up to (but not including) this time
        % step (k)
        for j = k-max(max(obj.Lags)):k-1
            col = j - (k-max(max(obj.Lags)))+1;
            for typidx = 1:size(obj.Types,2)
               for i = 1:size(obj.Measurements{typidx},2)
                   varidx = obj.Measurements{typidx}(i);
                   
                   % if measurements are in past data (ie real plant data)
                   % take it from there
                   if(j<=past_length)
                        measured_values{typidx}(varidx, col) = ...
                            past_data{typidx}(varidx,j);
                        
                   % if measurements needed at this time step were
                   % predicted at past time steps use these predictions (ie
                   % outputs predicted by the model)
                   elseif(ismember(varidx, obj.Predicted{typidx}))
                        measured_values{typidx}(varidx, col) = ...
                            predicted_values{typidx}(varidx, j-past_length);
                        
                   % Otherwise, data must have been given as future
                   % measurements (ie control input trajectory or output
                   % from another model)
                   else
                        measured_values{typidx}(varidx, col) = ...
                            future_data{typidx}(varidx, j-past_length);
                   end

               end
            end
        end
        
        %%%% use 'Step_Forward' to predict the values one step ahead %%%%
        
        % if this point in the code has been reached and the model is not
        % current the user has already selected to ignore the uncurent
        % model.  To avoid showing them this warning message repeatedly,
        % set current to true for now and then set it back laster.
        model_was_current = obj.Model_Current; 
        obj.Model_Current = true;
        try % so that model status will be properly updated even if an error occurs
            [step_results] = Step_Forward(obj, 'struct', measured_values);
        catch exception
            obj.Model_Current = model_was_current;
            rethrow(exception)
        end
        obj.Model_Current = model_was_current;
        
        %%%% store the resultant value as a predicted value for use in the
        %%%% next time step %%%%
        for typidx = 1:size(obj.Types,2)
           for i = 1:size(obj.Predicted{typidx},2)
               varidx = obj.Predicted{typidx}(i);
               
               predicted_values{typidx}(varidx, k-past_length) = ...
                   step_results{typidx}(varidx,1);
           end
        end
    end
end