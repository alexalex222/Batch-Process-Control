function [] = Fit_Model(obj, varargin)
   % Determine if necessairy model parameters have been set to allow
   % fitting
    if( ~obj.Measurements_set )
        error('Set measured variables using the ''AddMeasured'' method');
    elseif(~obj.Predicted_set)
        error('Set predicted variables using the ''AddPredicted'' method');
    elseif(~obj.Num_batch_set && obj.Is_Batch)
        error('Set number of batches using the ''SetBatchNumber'' method');
    elseif(~obj.Batch_time_set && obj.Is_Batch )
        error('Set batch time using the ''SetBatchLength'' method');
    end
    

   %%%% initialize struct to hold data %%%%
   raw_data = cell(size(obj.Types,2),2);  % Columns represent (1) fit or (2) validation
   
   
   %%%% Get data that was passed in %%%%
   if(mod(size(varargin,2),2) ~=0)
       error('There is a problem with the passed variables');
   end
   nvarpass = size(varargin,2)/2;
   
   for n = 1:nvarpass
       name = varargin{(n-1)*2+1};
       k = strfind(name,'_');
       
       % make sure variable is specified as either fit or validation 
       if(isempty(k) ||...
           ~(strcmp(name(k(end)+1:end),'fit') ||...
           strcmp(name(k(end)+1:end), 'val')))
       
           error(['''' name ''' must specify fit (''*_fit'') or validation (''*_val'')'])
       end
       type = name(1:k(end)-1);
       typidx = getTypeIndex(obj, type);
       
       % Determine if this data is fit or validation data
       if(strcmp(name(k(end)+1:end),'fit'))
           col = 1;  % i.e. Data is fit data
       else
           col = 2;  % i.e. Data is validation data
       end
       
       % store the data
       raw_data{typidx,col} = varargin{n*2};
   end
   
   % check to make sure all needed data was passed
   for typidx = 1:size(obj.Types,2)
      if(isempty(raw_data{typidx,1})&&...
              ~(isempty(obj.Measurements{typidx}) && isempty(obj.Predicted{typidx}))) 
         error(['The values for ''' obj.Types{typidx} '_fit'' data were not specified']);
      end
      if(isempty(raw_data{typidx,2})&&...
              ~(isempty(obj.Measurements{typidx}) && isempty(obj.Predicted{typidx}))) 
         error(['The values for ''' obj.Types{typidx} '_val'' data were not specified']);
      end
   end
   
   %%%% Confirm that the dimensions of the data added are correct %%%%
   for typidx = 1:size(obj.Types,2)
       if( ~(isempty(obj.Measurements{typidx}) && isempty(obj.Predicted{typidx})))
           if(size(raw_data{typidx,1},2) ~= obj.nFE)
              error(['The number of columns of ''' obj.Types{typidx}...
                  '_fit'' must be equal to the number of finite elements']);
           elseif(size(raw_data{typidx,2},2) ~= obj.nFE)
              error(['The number of columns of ''' obj.Types{typidx}...
                  '_val'' must be equal to the number of finite elements']);
           elseif(size(raw_data{typidx,1},3) ~= obj.nb_fit)
               error(['The size of the third dimension of ''' obj.Types{typidx}...
                   '_fit'' must be equal to the number of fit batches']);
           elseif(size(raw_data{typidx,2},3) ~= obj.nb_val)
               error(['The size of the third dimension of ''' obj.Types{typidx}...
                   '_val'' must be equal to the number of validation batches']);
           elseif(~isempty(obj.Measurements{typidx}) && ...
                   size(raw_data{typidx,1},1) < max(obj.Measurements{typidx}))
              error(['''' obj.Types{typidx} '_fit'' must have at least '...
                  num2str(max(obj.Measurements{typidx})) ' rows']); 
           elseif(~isempty(obj.Measurements{typidx}) && ...
                   size(raw_data{typidx,2},1) < max(obj.Measurements{typidx}))
              error(['''' obj.Types{typidx} '_val'' must have at least '...
                  num2str(max(obj.Measurements{typidx})) ' rows']); 
           elseif(~isempty(obj.Predicted{typidx}) && ...
                   size(raw_data{typidx,2},1) < max(obj.Predicted{typidx}))
               error(['''' obj.Types{typidx} '_val'' must have at least '...
                  num2str(max(obj.Predicted{typidx})) ' rows']); 
           elseif(~isempty(obj.Predicted{typidx}) && ...
                   size(raw_data{typidx,1},1) < max(obj.Predicted{typidx}))
               error(['''' obj.Types{typidx} '_fit'' must have at least '...
                  num2str(max(obj.Predicted{typidx})) ' rows']); 
           end
       end
   end
   
   %%%% Preallocate space for final data structure %%%%
   tot_measured =sum(sum(obj.Lags));
   tot_predict = 0;
   for typidx = 1:size(obj.Types,2)
        tot_predict = tot_predict + size(obj.Predicted{typidx},2);
   end
   
   %   notice that the maximum lag is subtracted because the model is not
   %   applicable until enough time steps have passed for the data to be
   %   lagged propperly    V  
   row_per_batch = obj.nFE-max(max(obj.Lags));
   tot_rows = row_per_batch*obj.nb_fit;
   
   indexBatchOne = max(max(obj.Lags))+1:1:obj.nFE;
   
   % rows denote times (stacked batches), columns dentoe variables
   measured_fit = zeros(obj.nb_fit*(row_per_batch),tot_measured); 
   predict_fit = zeros(size(measured_fit,1), tot_predict);
   
   measured_val = zeros(obj.nb_val*(row_per_batch),tot_measured);
   predict_val = zeros(size(measured_val,1), tot_predict);
   
   
   %%%% Unfold and lag data %%%%
   
   % Measurement fit data
   col = 1;
   row = 1;
   for typidx = 1:size(obj.Types,2)
       for i = 1:size(obj.Measurements{typidx},2)
           varidx = obj.Measurements{typidx}(i);
           for lag = 1:obj.Lags(typidx,varidx)
               for batch = 1:obj.nb_fit
                    measured_fit(row:row+row_per_batch-1,col) =...
                        raw_data{typidx,1}(varidx, indexBatchOne-lag, batch);
                    row = row+row_per_batch;
               end
               row = 1;
               col = col +1;
           end
       end
   end
   
   % Predicted fit data (no lags in predicted data)
   row =1;
   col = 1;
   for typidx = 1:size(obj.Types,2)
      for i = 1:size(obj.Predicted{typidx},2)
         varidx = obj.Predicted{typidx}(i);
         for batch = 1:obj.nb_fit
             predict_fit(row:row + row_per_batch-1, col) =...
                 raw_data{typidx,1}(varidx, indexBatchOne, batch);
             row = row+row_per_batch;
         end
         row = 1;
         col = col +1;
      end
   end
   
   
   %%%% Center and scale the data %%%%
   [measured_fit_scld, obj.Measure_Means, obj.Measure_Std] = zscore(measured_fit);
   [predict_fit_scld, obj.Predicted_Means, obj.Predicted_Std] = zscore(predict_fit);
   
   fprintf('> Fitting Model...\n');
   
   mean_rmse = cell(size(obj.Types,2), 1);
   mean_rmse_fit = cell(size(obj.Types,2), 1);
   coeffs = cell(numel(obj.Clusters));
   Centers = cell(numel(obj.Clusters));
   for i=1:numel(obj.Clusters)
        disp(['>> Progress - cluster ' num2str(obj.Clusters(i)) ' of ' num2str( obj.Clusters(end) )]);
        
        [mf, Centers{i}, dummy1, dummy2] =...
            Mult_Mod.fcm_wrapper(measured_fit_scld,obj.Clusters(i),obj.Cluster_Iter); 
        
        % Make weighted regressor matrix required for the regression:
        Xscld_T_wgtd = reshape(repmat(mf,1,(tot_measured+1))',tot_rows,(tot_measured+1)*obj.Clusters(i))...
            .*repmat([ones(tot_rows,1) measured_fit_scld],1,obj.Clusters(i));
        
        % IF using PLS, loop over each principal component possibility
        if(strcmpi(obj.Regression_Type, 'PLS'))
            nPC=obj.PLS_nPC;
        else
            nPC =1; % if not using PLS just fit the ols model once for each number of clusters
        end
        
        for k =1:numel(nPC)
            if(strcmpi(obj.Regression_Type, 'PLS'))
                % Check to make sure not using more PCs than varaiables
                if(size(Xscld_T_wgtd,2)<obj.PLS_nPC(k))
                    fprintf(['--- Could not fit a PLS model with ' num2str(obj.PLS_nPC(k)) ' PCs for this clustering\n']);
                    max_error(i,k) = inf;
                    continue;
                else
                    fprintf(['--- Fitting a PLS model with ' num2str(obj.PLS_nPC(k)) ' PCs\n']);
                    
                    %apply scaling to variables
                    [Xscld_T_wgtd_rescld, mu, stdev] = zscore(Xscld_T_wgtd);
                    
                    
                    [wstar{i,k},yl{i,k},xs,ys,coeffs{i,k},P{i,k}] = nipals_pls(Xscld_T_wgtd_rescld,predict_fit_scld, obj.PLS_nPC(k));
                    %[dummy1, dummy2, dummy3, dummy4, coeffs{i,k}] =  plsregress(Xscld_T_wgtd,predict_fit_scld, obj.PLS_nPC(k));
                    
                    beta_new = [mu./stdev*coeffs{i,k}; coeffs{i,k}./repmat(stdev',1,size(coeffs{i,k},2))];
                    
                    xhat = xs*P{i,k}';
                    err = Xscld_T_wgtd - xhat;
                    SPE = sum(err.^2,2);
                    
                    m = mean(SPE);
                    v = var(SPE);
                    
                    h(i,k) = 2*m^2/v;
                    g(i,k) = v/(2*m);
                    
                    numobs(i,k) = size(Xscld_T_wgtd,1);
                    
                    sa{i,k} = std(xs);
                end
            else
                % Linear regression and associated covariance of parameters:
                coeffs{i,k} = (Xscld_T_wgtd'*Xscld_T_wgtd)\Xscld_T_wgtd'*predict_fit_scld;
            end

            % *** Temporarily store the results for this number of clusters in
            % the model parameters so that other model methods can be used to
            % determine RMSE etc.  
            obj.Centers = Centers{i};
            obj.Opt_Clusters = obj.Clusters(i);
            obj.Coeffs = coeffs{i,k};
            obj.Model_Current = true;
            obj.Model_Fitted = true;

            % Try catch block so that if an error occurs all model status
            % parameters can be fixed (see above comment)
            try
                rmse_all = cell(size(obj.Types,2),obj.nb_val);
                for n =1:obj.nb_val
                    val_batch = cell(size(obj.Types,2),1);
                    for typidx = 1:size(obj.Types,2)
                        % Add measurement data
                        for j = 1:size(obj.Measurements{typidx},2) 
                            varidx = obj.Measurements{typidx}(j);
                            val_batch{typidx,1}(varidx,:) = raw_data{typidx,2}(varidx,:,n);
                        end
                        for j = 1:size(obj.Predicted{typidx},2)
                            varidx = obj.Predicted{typidx}(j);
                            % if data is not already added in the measurement
                            % data, add observed outputs so that error can be
                            % calculated
                            if(~ismember(varidx,obj.Measurements{typidx}))
                               val_batch{typidx,1}(varidx,:) = raw_data{typidx,2}(varidx,:,n); 
                            end

                        end
                    end
                    rmse_all(:,n) = obj.RMSE('struct', val_batch); 
                end
                
                rmse_fit_all = cell(size(obj.Types,2),obj.nb_fit);
                for n =1:obj.nb_fit
                    fit_batch = cell(size(obj.Types,2),1);
                    for typidx = 1:size(obj.Types,2)
                        % Add measurement data
                        for j = 1:size(obj.Measurements{typidx},2) 
                            varidx = obj.Measurements{typidx}(j);
                            fit_batch{typidx,1}(varidx,:) = raw_data{typidx,1}(varidx,:,n);
                        end
                        for j = 1:size(obj.Predicted{typidx},2)
                            varidx = obj.Predicted{typidx}(j);
                            % if data is not already added in the measurement
                            % data, add observed outputs so that error can be
                            % calculated
                            if(~ismember(varidx,obj.Measurements{typidx}))
                               fit_batch{typidx,1}(varidx,:) = fit_batch{typidx,1}(varidx,:,n); 
                            end

                        end
                    end
                    rmse_fit_all(:,n) = obj.RMSE('struct', fit_batch); 
                end


                for typidx = 1:size(obj.Types,2)
                    for j = 1:size(obj.Predicted{typidx},2)
                        varidx = obj.Predicted{typidx}(j);
                        tmp = cell2mat(rmse_all(typidx,:));
                        mean_rmse{typidx}(varidx,i,k) = mean(tmp(varidx,:));
                    end
                end
                
                for typidx = 1:size(obj.Types,2)
                    for j = 1:size(obj.Predicted{typidx},2)
                        varidx = obj.Predicted{typidx}(j);
                        tmp = cell2mat(rmse_fit_all(typidx,:));
                        mean_rmse_fit{typidx}(varidx,i,k) = mean(tmp(varidx,:));
                    end
                end


            catch exception
                obj.Centers = [];
                obj.Opt_Clusters = [];
                obj.Coeffs = [];
                obj.Model_Current = false;
                obj.Model_Fitted = false;
                rethrow(exception);
            end
            
            obj.Centers = [];
            obj.Opt_Clusters = [];
            obj.Coeffs = [];
            obj.Model_Current = false;
        
            tmp = cell2mat(mean_rmse);
            max_error(i,k) = max(tmp(:,i,k));  % find the maximum rmse (out of outputs predicted) for each number of clusters and principal components
        end 
   end
   
   
   [dummy opt_nPC_idx] = min(min(max_error));
   [dummy opt_clust_idx] = min(max_error(:,opt_nPC_idx));
   
   
   % PLS related results
   if(strcmpi(obj.Regression_Type, 'PLS'))
       obj.g = g(opt_clust_idx, opt_nPC_idx);
       obj.h = h(opt_clust_idx, opt_nPC_idx);
       obj.wstar = wstar{opt_clust_idx, opt_nPC_idx};
       obj.yl = yl{opt_clust_idx, opt_nPC_idx};
       obj.loading = P{opt_clust_idx, opt_nPC_idx};
       obj.numobs = numobs(opt_clust_idx, opt_nPC_idx);
       obj.sa = sa{opt_clust_idx, opt_nPC_idx};
   end
   
   obj.Opt_Clusters = obj.Clusters(opt_clust_idx);
   obj.Opt_nPC= obj.PLS_nPC(opt_nPC_idx);
   obj.Centers = Centers{opt_clust_idx};
   obj.Coeffs = coeffs{opt_clust_idx,opt_nPC_idx};
   for typidx = 1:size(obj.Types,2)
       for j = 1:size(obj.Predicted{typidx},2)
           varidx = obj.Predicted{typidx}(j);
          % obj.RMSE_val{typidx}(varidx,1) = mean_rmse{typidx}(varidx,opt_clust_idx, opt_nPC_idx);
          obj.RMSE_val{typidx}(varidx,1,:,:) = mean_rmse{typidx}(varidx,:, :);
       end
   end
   
   for typidx = 1:size(obj.Types,2)
       for j = 1:size(obj.Predicted{typidx},2)
           varidx = obj.Predicted{typidx}(j);
          % obj.RMSE_fit{typidx}(varidx,1) = mean_rmse_fit{typidx}(varidx,opt_clust_idx, opt_nPC_idx);
          obj.RMSE_fit{typidx}(varidx,1,:,:) = mean_rmse_fit{typidx}(varidx,:, :);
       end
   end
    
   obj.Model_Current = true;
   obj.Model_Fitted = true;
end
