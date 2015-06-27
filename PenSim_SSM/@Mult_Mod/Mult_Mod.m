classdef Mult_Mod < handle
   properties(GetAccess = 'public', SetAccess = 'private')
       %%%% Model Structure %%%%
       Measurements = {};  % Inputs to the model 
       Predicted = {}; % Variables that this model predicts (struct, index denotes type)
       Lags= []; % Lag Structure (array, row denotes type, column denotes index, value denotes lags)
       Clusters = 1:10; % Different numbers of clusters to try (ie 1:10 for between one and ten clusters)
       Is_Batch;  %Boolean.  Determines whether data is unfolded or not. 
       Regression_Type = 'OLS';   % Set to 'OLS' for ordinary Least Squares regression or PLS for partial least squares regression
       PLS_nPC = 4;
       
       %%%% Fitting Options %%%%
       Cluster_Iter = 15;
       
       %%%% Batch Information (only used if is batch) %%%% 
       nb_fit = -1;
       nb_val = -1;
       nFE= -1;
       samt = -1;
       
       
       %%%% Model Status %%%%
       Model_Current = false; % Set to true when the model has been fit to data and is ready to predict
                      %  (reset to fasle when any changes in the model
                      %  parameters are made)
       Model_Fitted = false; % set to true at the end of the 'Fit_Model' method               
       
       
       % Model parameters status
       Measurements_set = false;
       Predicted_set = false;
       Num_batch_set = false;
       Batch_time_set = false;
       
       %%%% Fit Results %%%%
       Measure_Means;
       Measure_Std;
       Predicted_Means;
       Predicted_Std;
       
       Opt_Clusters;
       Opt_nPC;
       Centers; % Clustering centers
       Coeffs;
       RMSE_val;
       RMSE_fit;
       
       % PLS Model
       loading = [];
       wstar = [];
       yl = [];
       g = [];
       h = [];
       numobs = [];
       sa = [];
       
   end
   
   properties(Constant)
      Types = {'state', 'output', 'input', 'other'}; 
   end
   
   properties(Dependent = true)
      NumMeasured
      NumPredicted
   end
   
   methods
      function value = get.NumMeasured(obj)
          value =0;
          for typidx = 1:size(obj.Types,2)
            value = size(obj.Measurements{typidx},2) + value;
          end
      end 
      
      function value = get.NumPredicted(obj)
          value =0;
          for typidx = 1:size(obj.Types,2)
            value = size(obj.Predicted{typidx},2) + value;
          end
          
      end
   end
   
   methods(Access = 'public')
       %%%% CONSTRUCTOR -> initialize variables and set flags %%%%
       function obj = Mult_Mod(Is_Batch)
           if(nargin < 1)
               Is_Batch = true;
           end
           
           obj.Is_Batch = Is_Batch;
           obj.Measurements = cell(size(obj.Types,2),1);
           obj.Predicted = cell(size(obj.Types,2),1);
       end
       
       %%%% Variable selection methods %%%%
       function [] = AddMeasured(obj, type, index, lag)
           typidx = getTypeIndex(obj, type);
           
           
           if (~ismember(index, obj.Measurements{typidx}))
               obj.Measurements{typidx} = [obj.Measurements{typidx} index];
               obj.Measurements{typidx} = sort(obj.Measurements{typidx});
           end
           obj.Lags(typidx, index) = lag;
           
           obj.Model_Current = false;
           obj.Measurements_set = true;
       end
       function [] = AddPredicted(obj, type, index)
          typidx = getTypeIndex(obj, type);
          
          if(~ismember(index, obj.Predicted{typidx}))
              obj.Predicted{typidx} = [obj.Predicted{typidx} index];
              obj.Predicted{typidx} = sort(obj.Predicted{typidx});
          end
          
          obj.Model_Current = false;
          obj.Predicted_set = true;
       end
       function [] = RemoveMeasured(obj, type, index)
          typidx = getTypeIndex(obj, type);
          
          % if the index is in the measurement set remove it and reset lags
          if(ismember(index, obj.Measurements{typidx}))
              obj.Measurements{typidx} = ...
                  obj.Measurements{typidx}(obj.Measurements{typidx}~=index);
              
              % reset lags
              obj.Lags(typidx,index) = 0;
          end
          
          obj.Model_Current = false;
       end
       function [] = RemovePredicted(obj,type, index)
          typidx = getTypeIndex(obj, type);
          
          if(ismember(index, obj.Predicted{typidx}))
             obj.Predicted{typidx} =...
                 obj.Predicted{typidx}(obj.Predicted{typidx}~=index);
              
          end
          
          obj.Model_Current = false;
       end
       
       
       
       %%%% Model Structure Parameters %%%%
       function [] = SetLag(obj, type, index, lag)
          typidx = getTypeIndex(obj, type);
          
          obj.Lags(typidx, index) = lag;
          
          obj.Model_Current = false;
       end
       function [] = SetClusters(obj, Clusters)
           obj.Clusters = Clusters;
           
           obj.Model_Current = false;
       end
       function [] = SetRegressionType(obj, type)
          obj.Regression_Type = type;
          
          obj.Model_Current = false;
       end
       function[] = SetPLSnPC(obj, nPC)
          obj.PLS_nPC = nPC;
          
          obj.Model_Current = false;
       end
       function [] = SetBatchLength(obj, num_FE, samp_time)
           obj.nFE = num_FE;
           obj.samt = samp_time;
           
           obj.Model_Current = false;
           obj.Batch_time_set = true;
       end
       function [] = SetBatchNumber(obj, nb_fit, nb_val)
           obj.nb_fit = nb_fit;
           obj.nb_val = nb_val;
           
           obj.Model_Current = false;
           obj.Num_batch_set = true;
       end
       function index = getTypeIndex(obj, type)
            index = find(strcmp(obj.Types,type));
            if(isempty(index))
                error(['The type "' type '" doesn''t exist']) 
            end
       end
       function [] = SetClusterIteration(obj, Cluster_Iter)
           obj.Cluster_Iter = Cluster_Iter;
           
           obj.Model_Current = false;
       end
       
       %%%% Model Fitting Methods %%%%
       [] = Fit_Model(obj, varargin);
       
       
       RMSE = RMSE(obj, varargin);
       
       %%%% model validity limits for PLS models
       SPE_lim = SPE_lim(obj, Prob);
       Tsqrd_lim = Tsqrd_lim(obj, Prob);
       [SPE Tsqrd] = getSPEandTsqrd(obj, varargin);
       
       %%%% Model Application %%%%
       Step_Forward = Step_Forward(obj, varargin);
       Predict_Forward = Predict_Forward(obj, past_length, steps, varargin);
       
       
       %%%% Display functions %%%%
       Plot_Traj = Plot_Trag(obj,varargin);
       Disp_RMSE_val = Disp_RMSE_val(obj)
       
       function [] = Disp_Struct(obj)
          fprintf('------------------------------------------------------\n');
          
          
          fprintf('Clusters:\n');
          fprintf(['   [' num2str(obj.Clusters) ']\n\n']);
           
           fprintf('Measured Variables:\n');
          for i = 1:size(obj.Types,2)
             if(~isempty(obj.Measurements{i}))
                for j = 1:size(obj.Measurements{i},2)
                    fprintf(['   ' obj.Types{i} ' ' num2str(obj.Measurements{i}(j))...
                        ', \tL = ' num2str(obj.Lags(i,obj.Measurements{i}(j))) '\n']);
                end
             end
          end
          fprintf('\nPredicted Variables:\n')
          for i = 1:size(obj.Types,2)
              if(~isempty(obj.Predicted{i}))
                  for j = 1:size(obj.Predicted{i},2)
                     fprintf(['   ' obj.Types{i} ' ' num2str(obj.Predicted{i}(j)) '\n']); 
                  end
              end
          end
       end
   end 
   
   methods(Static = true, Access = 'private')
       [MF C centers J] = fcm_wrapper(x,L,nFMCiter);
       [center, U, obj_fcn] = fcm(data, cluster_n, options);
       U = initfcm(cluster_n, data_n);
       [U_new, center, obj_fcn] = stepfcm(data, U, cluster_n, expo);
       out = distfcm(center, data);
   end
end

