clear variables;
close all;

% load 'bdb_v1_3.mat' 
load 'test9' 

fname = ['Q_model_SSM'];

% Principal components to try
PCs_to_try = 2:1:30;

% whether to include inputs in the model or not
use_input = true;

%set number of fitting, validation, testing
nb_fit = 40;
nb_val = 10;
nb_test =10;
nb_PRBS = 0;% bdb.nb_PRBS;

% Allow us to select the desired initial state mesasuremetns
desired_initial_states = [1 2 3 4 5 6 7];

% allow us to select which U and Y variabels to use
desired_outputs = [ 2 3 4 5]; % Select output variables we want to use
desired_inputs = []; % Select input variables we want to use
desired_qualities = [1 2 3]; % Select quality variables we want to predict

%get data from database and initial conditions

ynom = bdb.data.ynom(:,:);
unom = bdb.data.unom(:,:);
initialcond_nom = bdb.data.xnom(:,1);

%Create temporary matrix for y_fit, u_fit, q_fit and initialcond_fit
y_fit=zeros(size(bdb.data.y_PRBS,1),size(bdb.data.y_PRBS,2),nb_fit);
u_fit=zeros(size(bdb.data.u_PRBS,1),size(bdb.data.u_PRBS,2),nb_fit);
q_fit = zeros(size(bdb.data.q_PRBS,1),size(bdb.data.q_PRBS,2),nb_fit);
initialcond_fit = zeros(nb_fit,size(bdb.data.x_PRBS,1));

%Create temporary matrix for y_fit, u_fit, q_fit and initialcond_fit
y_fit=zeros(size(bdb.data.y_PRBS,1),size(bdb.data.y_PRBS,2),nb_fit);
u_fit=zeros(size(bdb.data.u_PRBS,1),size(bdb.data.u_PRBS,2),nb_fit);
q_fit = zeros(size(bdb.data.q_PRBS,1),size(bdb.data.q_PRBS,2),nb_fit);
initialcond_fit = zeros(nb_fit,size(bdb.data.x_PRBS,1));

%Get the 3D Matrix (variable x time x batch ) that will be used to fit the
%model
for i=1:nb_fit
    if(i<=nb_PRBS)
        y_fit(:,:,i) = bdb.data.y_PRBS(:,:,i);
        u_fit(:,:,i) = bdb.data.u_PRBS(:,:,i);
        q_fit(:,:,i) = bdb.data.q_PRBS(:,:,i);
        initialcond_fit(i,:) = bdb.data.x_PRBS(:,1,i);
    else
        y_fit(:,:,i) = bdb.data.y(:,:,i-nb_PRBS);
        u_fit(:,:,i) = bdb.data.u(:,:,i-nb_PRBS);
        q_fit(:,:,i) = bdb.data.q(:,:,i-nb_PRBS);
        initialcond_fit(i,:) = bdb.data.x(:,1,i-nb_PRBS);
    end
end

%Create temporary matrix for y_val, u_val, q_val and initialcond_val
y_val=zeros(size(bdb.data.y,1),size(bdb.data.y,2),nb_val);
u_val=zeros(size(bdb.data.u,1),size(bdb.data.u,2),nb_val);
q_val=zeros(size(bdb.data.q,1),size(bdb.data.q,2),nb_val);
initialcond_val = zeros(nb_val,size(bdb.data.x,1));

%Get the 3D Matrix (variable x time x batch ) that will be used to validate the
%model
for i=nb_fit+1:1:nb_val+nb_fit
    y_val(:,:,i-nb_fit) = bdb.data.y(:,:,i);
    u_val(:,:,i-nb_fit) = bdb.data.u(:,:,i);
    q_val(:,:,i-nb_fit) = bdb.data.q(:,:,i);
    initialcond_val(i-nb_fit,:) = bdb.data.x(:,1,i);
end


%Create temporary matrix for y_test, u_test, q_test and initialcond_test
y_test=zeros(size(bdb.data.y,1),size(bdb.data.y,2),nb_test);
u_test=zeros(size(bdb.data.u,1),size(bdb.data.u,2),nb_test);
q_test=zeros(size(bdb.data.q,1),size(bdb.data.q,2),nb_test);
initialcond_test = zeros(nb_test,size(bdb.data.x,1));

%Get the 3D Matrix (variable x time x batch ) that will be used to test the
%model
for i=(nb_fit+nb_val+1):1:nb_test+nb_fit+nb_val
    y_test(:,:,i-(nb_fit+nb_val)) = bdb.data.y(:,:,i);
    u_test(:,:,i-(nb_fit+nb_val)) = bdb.data.u(:,:,i);
    q_test(:,:,i-(nb_fit+nb_val)) = bdb.data.q(:,:,i);
    initialcond_test(i-(nb_fit+nb_val),:) = bdb.data.x(:,1,i);
end

%Select only the desired initial states
for i=1:numel(desired_initial_states)
    
initialcond_fit(:,i) = initialcond_fit(:,desired_initial_states(i));
initialcond_val(:,i) = initialcond_val(:,desired_initial_states(i));
initialcond_test(:,i) = initialcond_test(:,desired_initial_states(i));
initialcond_nom(i) = initialcond_nom(desired_initial_states(i));

end

%Delete the undesired initial conditions
for i=size(bdb.data.x,1):-1:numel(desired_initial_states)+1
    
initialcond_fit(:,i) = [];
initialcond_val(:,i) = [];
initialcond_test(:,i) = [];
initialcond_nom(i) = [];
    
end

%align input by adding extra column of input
u_fit(:,size(u_fit,2)+1,:) = u_fit(:,size(u_fit,2),:);
u_val(:,size(u_val,2)+1,:) = u_val(:,size(u_val,2),:);
u_test(:,size(u_test,2)+1,:) = u_test(:,size(u_test,2),:);

unom(:,size(unom,2)+1) = unom(:,size(unom,2));

% Create temporary matrix to allocate concat_fit, concat_val, concat_test
% and concat_nom
concat_fit = zeros(numel(desired_outputs)+numel(desired_inputs),size(y_fit,2),size(y_fit,3));
concat_val = zeros(numel(desired_outputs)+numel(desired_inputs),size(y_val,2),size(y_val,3));
concat_test = zeros(numel(desired_outputs)+numel(desired_inputs),size(y_test,2),size(y_test,3));
concat_nom = zeros(numel(desired_outputs)+numel(desired_inputs),size(ynom,2));
 
index = 1;
for i=1:numel(desired_outputs)
    
        concat_fit(index,:,:) = y_fit(desired_outputs(i),:,:);
        concat_val(index,:,:) = y_val(desired_outputs(i),:,:);
        concat_test(index,:,:) = y_test(desired_outputs(i),:,:);
        concat_nom(index,:) = ynom(desired_outputs(i),:);
        
   if i==numel(desired_outputs)
       
        if(use_input)
            
            for j=1:numel(desired_inputs)
            concat_fit(index+1,:,:) = u_fit(desired_inputs(j),:,:);
            concat_val(index+1,:,:) = u_val(desired_inputs(j),:,:);
            concat_test(index+1,:,:) = u_test(desired_inputs(j),:,:);
            concat_nom(index+1,:) = unom(desired_inputs(j),:);
            
            index = index + 1;
            
            end
            
        end
        
    end
    index = index + 1;
end

%Create temporary matrix to allocate xnom
xnom = zeros(1,size(concat_fit,2)*size(concat_fit,1));

%unfolding the 3D array to 2D matrices
index = 1;
for j=1:size(concat_fit,2)
    for i = 1:size(concat_fit,1)
        xnom(index)=concat_nom(i,j);
        index = index + 1;
    end
end

%Unfold the fit data in x_fit

x_fit=reshape(concat_fit, size(concat_fit,1)*size(concat_fit,2),size(concat_fit,3));
x_fit=x_fit';


% Unfold val data in x_val

x_val=reshape(concat_val,size(concat_val,1)*size(concat_val,2),size(concat_val,3));
x_val=x_val';


% Unfold test data in x_test

x_test=reshape(concat_test, size(concat_test,1)*size(concat_test,2),size(concat_test,3));
x_test=x_test';


% Unfold the quality data and select the actual available measurement
qual_fit=q_fit(desired_qualities,end,:);
qual_fit=reshape(qual_fit,size(q_fit,1),size(q_fit,3));
qual_fit=qual_fit';

qual_val=q_val(desired_qualities,end,:);
qual_val=reshape(qual_val,size(q_val,1),size(q_val,3));
qual_val=qual_val';

qual_test=q_test(desired_qualities,end,:);
qual_test=reshape(qual_test,size(q_test,1),size(q_test,3));
qual_test=qual_test';


%Put the initial states in the beginning of the matrix
x_fit = [initialcond_fit(:,1:size(initialcond_fit,2)) x_fit];
x_val = [initialcond_val(:,1:size(initialcond_val,2)) x_val];
x_test = [initialcond_test(:,1:size(initialcond_test,2)) x_test];
xnom = [initialcond_nom(1:size(initialcond_nom,1))' xnom];


[x_fit_centred, x_fit_mu, x_fit_sigma] = zscore(x_fit);

x_val_centred = (x_val - repmat(x_fit_mu,size(x_val,1), 1)) ./repmat(x_fit_sigma,size(x_val,1),1);
[qual_fit_centred, qual_fit_mu, qual_fit_sigma]  = zscore(qual_fit);
qual_val_centred = (qual_val - repmat(qual_fit_mu,size(qual_val,1), 1)) ./repmat(qual_fit_sigma,size(qual_val,1),1);
x_test_centred = (x_test - repmat(x_fit_mu,size(x_test,1), 1)) ./repmat(x_fit_sigma,size(x_test,1),1);
qual_test_centred = (qual_test - repmat(qual_fit_mu,size(qual_test,1), 1)) ./repmat(qual_fit_sigma,size(qual_test,1),1);

% UP TO HERE IT IS ALL ABOUT UNFOLDING AND SCALING

%AFTER HERE WE CARRY OUT THE PLS AND SAVE THE VALUES


for i = 1:size(PCs_to_try,2);
    [xl,yl,xs,ys,beta, P] = nipals_pls(x_fit_centred,qual_fit_centred,PCs_to_try(i));
    beta_stored(:,:,i) = beta;
    yhat = x_val_centred*beta;
    error = (yhat(:,2:3)-qual_val_centred(:,2:3)).^2;
    sumerror(i) = sum(sum(error));
end

[wstarAll, ~, ~, ~, ~, PAll] = nipals_pls(x_fit_centred, qual_fit_centred, size(x_fit_centred,1));





