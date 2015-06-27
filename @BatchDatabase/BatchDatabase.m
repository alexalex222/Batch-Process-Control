classdef BatchDatabase < handle
   properties(GetAccess = public, SetAccess = private)
        nb; % Number of batches in the database
        nb_PRBS; % Number of PRBS batches in the database
        nx; % Number of states;
        ny; % Number of outputs;
        nu; % Number of inputs;     
        nq; % Number of quality measurements;
        
        samt; % Sampling time
        Tf;   % Final time
        
        data; % Batch data
       
        meas_std; % Standard devation in measurements
        
        ic_std; % Standard deviation in initial conditions 
   end
   
   properties(Dependent= true)
       nFE;
   end
   
   %------------------------------------------------------------------
   % Dependant Properties Methods
   methods
      function value = get.nFE(obj)
          value = obj.Tf/obj.samt;
       end 
   end
   
   methods(Access= public)
       % Constructor:
       function obj= BatchDatabase(nx, ny, nu, nq, Tf, samt)
           % Allow user to specify number of states, outputs, inputs, and qualities for their batch
          obj.nx = nx; 
          obj.ny = ny;
          obj.nu = nu;
          obj.nq = nq;
          obj.samt = samt;
          obj.Tf = Tf;
          
          % Initialize the database with no batches
          obj.nb = 0;
          obj.nb_PRBS = 0;
          
          % Initialize measurement standard devaition to none
          obj.meas_std = zeros(ny,1);
          
          % Initialize inital condition standard devation to none
          obj.ic_std = zeros(nx,1);
          
       end
       
       % Add measurement variance
       function obj = AddMeasStd(obj, meas_std)
           if(size(meas_std,1)~=obj.ny)
              error(['The number of rows in meas_var must be ny:' num2str(obj.ny)]);
           end
           if(size(meas_std,2) ~=1)
              error('meas_var must be a column vector');
           end
          obj.meas_std = meas_std;
       end
       
       % Add initial condition standard devation 
       function [] = AddICStd(obj, ic_std)
           if(size(ic_std,1)~=obj.nx)
              error(['The number of rows in ic_std must be nx:' num2str(obj.nx)]);
           end
           if(size(ic_std,2) ~=1)
              error('ic_std must be a column vector');
           end
           obj.ic_std = ic_std;
       end
       
       % Add nominal batch
       function [] = AddNomBatch(obj, xnom, ynom, unom, qnom)
           % Throw errors for sampling time issues
          if(size(xnom,2)~= obj.nFE+1)
              error(['The number of rows in xnom must be ' num2str(obj.nFE +1)]);
          elseif(size(ynom,2)~= obj.nFE+1)
              error(['The number of rows in ynom must be ' num2str(obj.nFE +1)]);
          elseif(size(unom,2)~= obj.nFE)
              error(['The number of rows in unom must be ' num2str(obj.nFE)]);
          elseif(size(qnom,2)~= obj.nFE+1)
              error(['The number of rows in qnom must be ' num2str(obj.nFE +1)]);
          end
          
          % Throw errors for number of measurement issues
          if(size(xnom,1)~= obj.nx)
              error(['The number of columns in xnom must be ' num2str(obj.nx)]);
          elseif(size(ynom,1)~= obj.ny)
              error(['The number of columns in ynom must be ' num2str(obj.ny)]);
          elseif(size(unom,1)~= obj.nu)
              error(['The number of columns in unom must be ' num2str(obj.nu)]);
          elseif(size(qnom,1)~= obj.nq)
              error(['The number of columns in qnom must be ' num2str(obj.nq)]);
          end
           
          obj.data.xnom(:,:) = xnom;
          obj.data.ynom(:,:) = ynom;
          obj.data.unom(:,:) = unom;
          obj.data.qnom(:,:) = qnom;
       end
       
       
       % Add a new batch observation to the database
       function [] = AddBatch(obj, x, y, u, q)
           % Throw errors for sampling time issues
          if(size(x,2)~= obj.nFE+1)
              error(['The number of rows in x must be ' num2str(obj.nFE +1)]);
          elseif(size(y,2)~= obj.nFE+1)
              error(['The number of rows in y must be ' num2str(obj.nFE +1)]);
          elseif(size(u,2)~= obj.nFE)
              error(['The number of rows in u must be ' num2str(obj.nFE)]);
          elseif(size(q,2)~= obj.nFE+1)
              error(['The number of rows in q must be ' num2str(obj.nFE +1)]);
          end
          
          % Throw errors for number of measurement issues
          if(size(x,1)~= obj.nx)
              error(['The number of columns in x must be ' num2str(obj.nx)]);
          elseif(size(y,1)~= obj.ny)
              error(['The number of columns in y must be ' num2str(obj.ny)]);
          elseif(size(u,1)~= obj.nu)
              error(['The number of columns in u must be ' num2str(obj.nu)]);
          elseif(size(q,1)~= obj.nq)
              error(['The number of columns in q must be ' num2str(obj.nq)]);
          end
           
          obj.data.x(:,:, obj.nb+1) = x;
          obj.data.y(:,:, obj.nb+1) = y;
          obj.data.u(:,:, obj.nb+1) = u;
          obj.data.q(:,:, obj.nb+1) = q;
          
          obj.nb = obj.nb+1; 
       end
       
       % Add PRBS batch to database
       function [] = AddPRBSBatch(obj, x, y, u, q, Tsw)  
           % Tsw = average switching time
           
           % Throw errors for sampling time issues
          if(size(x,2)~= obj.nFE+1)
              error(['The number of rows in x must be ' num2str(obj.nFE +1)]);
          elseif(size(y,2)~= obj.nFE+1)
              error(['The number of rows in y must be ' num2str(obj.nFE +1)]);
          elseif(size(u,2)~= obj.nFE)
              error(['The number of rows in u must be ' num2str(obj.nFE)]);
          elseif(size(q,2)~= obj.nFE+1)
              error(['The number of rows in q must be ' num2str(obj.nFE +1)]);
          end
          
          % Throw errors for number of measurement issues
          if(size(x,1)~= obj.nx)
              error(['The number of columns in x must be ' num2str(obj.nx)]);
          elseif(size(y,1)~= obj.ny)
              error(['The number of columns in y must be ' num2str(obj.ny)]);
          elseif(size(u,1)~= obj.nu)
              error(['The number of columns in u must be ' num2str(obj.nu)]);
          elseif(size(q,1)~= obj.nq)
              error(['The number of columns in q must be ' num2str(obj.nq)]);
          end
           
          obj.data.x_PRBS(:,:, obj.nb_PRBS +1) = x;
          obj.data.y_PRBS(:,:, obj.nb_PRBS +1) = y;
          obj.data.u_PRBS(:,:, obj.nb_PRBS +1) = u;
          obj.data.q_PRBS(:,:, obj.nb_PRBS +1) = q;
          obj.data.Tsw(obj.nb_PRBS +1) = Tsw;
          
          obj.nb_PRBS = obj.nb_PRBS +1;
       end
       
       % Genereate noise trajectory
       function noise = GenNoise(obj, seed)
          randn('seed', seed);
          
          noise= randn(obj.ny, obj.nFE+1);
   
          % Adjust to correct standard deviation (based on the user determined sensor variance)
          noise = noise.*repmat(obj.meas_std,1,obj.nFE+1); 
       end
       
       % Returns a matrix of inital conditions 
       function initial_conds = GetICs(obj, num_conds, seed)
           rand('seed',seed);
           
           initial_conds = 2.*rand(num_conds, obj.nx)-1;
           initial_conds = repmat(obj.data.xnom(:,1)', num_conds,1) +initial_conds.*repmat(obj.ic_std', num_conds,1);
       end
       
       
       % Plot trajectories
       function [] = PlotTraj(obj, type, index, batches, linespec)
           hold_state_orig = ishold();
           switch(type)
               case 'state'
                   to_plot = obj.data.x(index,:,batches);
                    
               case 'output'
                   to_plot = obj.data.y(index,:,batches);
                   
               case 'input'
                   to_plot = obj.data.u(index,:,batches);
                   
               case 'state_PRBS'
                   to_plot = obj.data.x_PRBS(index,:,batches);
                    
               case 'output_PRBS'
                   to_plot = obj.data.y_PRBS(index,:,batches);
                   
               case 'input_PRBS'
                   to_plot = obj.data.u_PRBS(index,:,batches);
           end
           
           gca;
           hold on;
           if strcmp(type, 'input') || strcmp(type, 'input_PRBS')
               for i = 1:size(to_plot,3)
                  stairs(0:obj.nFE-1, to_plot(:,:,i), linespec);
               end
           else
               for i = 1:size(to_plot,3)
                  plot(0:obj.nFE, to_plot(:,:,i), linespec);
               end
           end
           
           if(hold_state_orig)
              hold on;
           else
              hold off;
           end
       end
   end
end