clear variables;
close all;


load 'test4' 

fname = ['Q_model_final'];

model_use = 'MultiMod'; % either 'MissingData' or 'MultiMod'

use_input = true;  % whether to include inputs in the model or not

% Principal components to try
PCs_to_try = 2:1:30;

%set number of fitting, validation, testing
nb_fit = 40;
nb_val = 10;
nb_test =10;
nb_PRBS = 0;% bdb.nb_PRBS;


% Allow us to select the desired initial state mesasuremetns
% desired_initial_states = [1 2 3 4 5 6 7];
desired_initial_states = [1 3 4];

% allow us to select which U and Y variabels to use
desired_outputs = [ 2 3 4 5]; % Select output variables we want to use
desired_inputs = []; % Select input variables we want to use

desired_qualities = [1 2 3]; % Select quality variables we want to predict

% Select the index (time stamp) of available intermediate quality
IndexIni = 401;
TrajLength = 201;


%get data from database and initial conditions

ynom = bdb.data.ynom(:,:);
unom = bdb.data.unom(:,:);
initialcond_nom = bdb.data.xnom(:,1);

%Create temporary matrix for y_fit, u_fit, q_fit and initialcond_fit
% y_fit=zeros(size(bdb.data.y_PRBS,1),size(bdb.data.y_PRBS,2),nb_fit);
% u_fit=zeros(size(bdb.data.u_PRBS,1),size(bdb.data.u_PRBS,2),nb_fit);
% q_fit = zeros(size(bdb.data.q_PRBS,1),size(bdb.data.q_PRBS,2),nb_fit);

y_fit=zeros(size(bdb.data.y_PRBS,1),TrajLength,nb_fit);
u_fit=zeros(size(bdb.data.u_PRBS,1),TrajLength-1,nb_fit);
q_fit = zeros(size(bdb.data.q_PRBS,1),TrajLength,nb_fit);


initialcond_fit = zeros(nb_fit,size(bdb.data.x_PRBS,1));

%Get the 3D Matrix (variable x time x batch ) that will be used to fit the
%model
for i=1:nb_fit
    if(i<=nb_PRBS)
        y_fit(:,:,i) = bdb.data.y_PRBS(:,IndexIni:end,i);
        u_fit(:,:,i) = bdb.data.u_PRBS(:,IndexIni:end,i);
        q_fit(:,:,i) = bdb.data.q_PRBS(:,IndexIni:end,i);
        initialcond_fit(i,:) = bdb.data.x_PRBS(:,IndexIni,i);
    else
        y_fit(:,:,i) = bdb.data.y(:,IndexIni:end,i-nb_PRBS);
        u_fit(:,:,i) = bdb.data.u(:,IndexIni:end,i-nb_PRBS);
        q_fit(:,:,i) = bdb.data.q(:,IndexIni:end,i-nb_PRBS);
        initialcond_fit(i,:) = bdb.data.x(:,IndexIni,i-nb_PRBS);
    end
end

%Create temporary matrix for y_val, u_val, q_val and initialcond_val
% y_val=zeros(size(bdb.data.y,1),size(bdb.data.y,2),nb_val);
% u_val=zeros(size(bdb.data.u,1),size(bdb.data.u,2),nb_val);
% q_val=zeros(size(bdb.data.q,1),size(bdb.data.q,2),nb_val);

y_val=zeros(size(bdb.data.y,1),TrajLength,nb_val);
u_val=zeros(size(bdb.data.u,1),TrajLength-1,nb_val);
q_val=zeros(size(bdb.data.q,1),TrajLength,nb_val);

initialcond_val = zeros(nb_val,size(bdb.data.x,1));

%Get the 3D Matrix (variable x time x batch ) that will be used to validate the
%model
for i=nb_fit+1:1:nb_val+nb_fit
    y_val(:,:,i-nb_fit) = bdb.data.y(:,IndexIni:end,i);
    u_val(:,:,i-nb_fit) = bdb.data.u(:,IndexIni:end,i);
    q_val(:,:,i-nb_fit) = bdb.data.q(:,IndexIni:end,i);
    initialcond_val(i-nb_fit,:) = bdb.data.x(:,IndexIni,i);
end

%Create temporary matrix for y_test, u_test, q_test and initialcond_test
% y_test=zeros(size(bdb.data.y,1),size(bdb.data.y,2),nb_test);
% u_test=zeros(size(bdb.data.u,1),size(bdb.data.u,2),nb_test);
% q_test=zeros(size(bdb.data.q,1),size(bdb.data.q,2),nb_test);

y_test=zeros(size(bdb.data.y,1),TrajLength,nb_test);
u_test=zeros(size(bdb.data.u,1),TrajLength-1,nb_test);
q_test=zeros(size(bdb.data.q,1),TrajLength,nb_test);

initialcond_test = zeros(nb_test,size(bdb.data.x,1));

%Get the 3D Matrix (variable x time x batch ) that will be used to test the
%model
for i=(nb_fit+nb_val+1):1:nb_test+nb_fit+nb_val
    y_test(:,:,i-(nb_fit+nb_val)) = bdb.data.y(:,IndexIni:end,i);
    u_test(:,:,i-(nb_fit+nb_val)) = bdb.data.u(:,IndexIni:end,i);
    q_test(:,:,i-(nb_fit+nb_val)) = bdb.data.q(:,IndexIni:end,i);
    initialcond_test(i-(nb_fit+nb_val),:) = bdb.data.x(:,IndexIni,i);
end

% WHY DELETE THESE COLUMNS ?
%set initial conditions of some states to zero
% initialcond_fit(:,3:1:8) = [];
% initialcond_val(:,3:1:8) = [];
% initialcond_test(:,3:1:8) = [];
% initialcond_nom(3:1:8) = [];

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
% u_fit(:,241,:) = u_fit(:,240,:);
% u_val(:,241,:) = u_val(:,240,:);
% u_test(:,241,:) = u_test(:,240,:);
% 
% unom(:,241) = unom(:,240);

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



% Normalize the data
switch model_use
    case 'MultiMod'
        [x_fit_centred, x_fit_mu, x_fit_sigma] = zscore(x_fit);
    case 'MissingData'
%         x_fit_mu = xnom;
%         x_fit_sigma = std(x_fit - repmat(x_fit_mu, size(x_fit,1), 1 ));
%         x_fit_sigma( x_fit_sigma==0)=1;
%         x_fit_centred = (x_fit - repmat(x_fit_mu,size(x_fit,1), 1)) ./repmat(x_fit_sigma,size(x_fit,1),1);
        
        [x_fit_centred, x_fit_mu, x_fit_sigma] = zscore(x_fit);
        
        x_fit_sigma(x_fit_sigma==0) = 1;
end

x_val_centred = (x_val - repmat(x_fit_mu,size(x_val,1), 1)) ./repmat(x_fit_sigma,size(x_val,1),1);
[qual_fit_centred, qual_fit_mu, qual_fit_sigma]  = zscore(qual_fit);
qual_val_centred = (qual_val - repmat(qual_fit_mu,size(qual_val,1), 1)) ./repmat(qual_fit_sigma,size(qual_val,1),1);
x_test_centred = (x_test - repmat(x_fit_mu,size(x_test,1), 1)) ./repmat(x_fit_sigma,size(x_test,1),1);
qual_test_centred = (qual_test - repmat(qual_fit_mu,size(qual_test,1), 1)) ./repmat(qual_fit_sigma,size(qual_test,1),1);

%x_val_centred = [(ones(nb_val,1)) x_val_centred];
%x_test_centred = [(ones(nb_test,1)) x_test_centred];
if(strcmp('MultiMod', model_use))
    x_fit_centred(:,4:1:end) = x_fit_centred(:,4:1:end)./sqrt(size(x_fit_centred,2));
    x_val_centred(:,4:1:end) = x_val_centred(:,4:1:end)./sqrt(size(x_fit_centred,2));
    x_test_centred(:,4:1:end) = x_test_centred(:,4:1:end)./sqrt(size(x_fit_centred,2));
end


% UP TO HERE IT IS ALL ABOUT UNFOLDING AND SCALING

%AFTER HERE WE CARRY OUT THE PLS AND SAVE THE VALUES


for i = 1:size(PCs_to_try,2);
    [xl,yl,xs,ys,beta, P] = nipals_pls(x_fit_centred,qual_fit_centred,PCs_to_try(i));
    beta_stored(:,:,i) = beta;
    yhat = x_val_centred*beta;
    error = (yhat(:,2:3)-qual_val_centred(:,2:3)).^2;
    sumerror(i) = sum(sum(error));
end

[wstarAll, dummy1, dummy2, dummy3, dummy4, PAll] = nipals_pls(x_fit_centred, qual_fit_centred, size(x_fit_centred,1));


% % scl_mtx = ones(1,size(x_fit_centred,2));
% % 
% % for i =4:4:size(x_fit_centred,2)
% %    scl_mtx(i) = scl_mtx(i)*100; 
% % end
% 
% 
% x_fit_centred = x_fit_centred.*repmat(scl_mtx, nb_fit, 1);

[minsSPE ,nPC_idx] = min(sumerror);

[xl,yl,xs,ys,beta, P,dummy1,dummy2, W] = nipals_pls(x_fit_centred,qual_fit_centred,PCs_to_try(nPC_idx));
index = 1;
for i = 1:2
    for j = 1:size(concat_fit,1)
        index = 1;
        for k = j+(size(concat_fit,1))-1:size(concat_fit,1):(size(concat_fit,2)*size(concat_fit,1))+(j-1)
            xloading(index,j,i) = xl(k,i);
            index = index +1;
        end
    end
end

figure
subplot(2,4,1);
bar(xloading(:,1,1))
subplot(2,4,2);
bar(xloading(:,2,1))
subplot(2,4,3);
bar(xloading(:,3,1))
if(use_input)
    subplot(2,4,4);
    bar(xloading(:,4,1))
end
subplot(2,4,5);
bar(xloading(:,1,2))
subplot(2,4,6);
bar(xloading(:,2,2))
subplot(2,4,7);
bar(xloading(:,3,2))
if(use_input)
    subplot(2,4,8);
    bar(xloading(:,4,2))
end

figure
% subplot(1,2,1);
scatter(xs(:,1),xs(:,2))
% plot T^2 intervals
T_sqrd_95_2comp = (nb_fit-1)*(nb_fit+1)*2/(nb_fit*(nb_fit-2))*finv(.95,2,nb_fit-2);

s1 = std(xs(:,1));
s2 = std(xs(:,2));

xvals = -sqrt(T_sqrd_95_2comp)*s1:.01:sqrt(T_sqrd_95_2comp)*s1;
yvals = sqrt((T_sqrd_95_2comp-xvals.^2./s1^2)*s2^2);

hold on
plot(xvals, yvals, ':k');
yvals = -yvals;
plot(xvals, yvals, ':k');

T_sqrd_99_2comp = (nb_fit-1)*(nb_fit+1)*2/(nb_fit*(nb_fit-2))*finv(.99,2,nb_fit-2);

xvals = -sqrt(T_sqrd_99_2comp)*s1:.01:sqrt(T_sqrd_99_2comp)*s1;
yvals = sqrt((T_sqrd_99_2comp-xvals.^2./s1^2)*s2^2);

hold on
plot(xvals, yvals, '--k');
yvals = -yvals;
plot(xvals, yvals, '--k');


qual_pred_centred = x_test_centred*beta_stored(:,:,nPC_idx);

val_xscores = x_val_centred*xl;
fit_xscores = x_fit_centred*xl;



for i=1:size(qual_pred_centred,1)
    for j = 1:size(qual_pred_centred,2)
        qual_pred(i,j) = qual_pred_centred(i,j)*qual_fit_sigma(j)+qual_fit_mu(j);
        error(i,j) = abs(qual_pred(i,j) - qual_test(i,j));
    end
end
% RMSE = sqrt(sum(error.^2));

RMSE=zeros(1,3);
R2=zeros(1,3);
MAPE=zeros(1,3);

% for i=1:size(qual_pred_centred,2)
%     [R2(i),RMSE(i)]=rsquare(qual_test(:,i),qual_pred(:,i));
%     MAPE(i)=mean(abs(qual_test(:,i)-qual_pred(:,i))./qual_test(:,i));
% end



% fprintf(['RMSE:\n Glucose: ' num2str(RMSE(1)) '\n Biomass: ' num2str(RMSE(2)) '\n Penicillin: ' num2str(RMSE(3)) '\n']);
% fprintf(['MAPE:\n Glucose: ' num2str(MAPE(1)) '\n Biomass: ' num2str(MAPE(2)) '\n Penicillin: ' num2str(MAPE(3)) '\n']);
% fprintf(['R2:\n Glucose: ' num2str(R2(1)) '\n Biomass: ' num2str(R2(2)) '\n Penicillin: ' num2str(R2(3)) '\n']);

%%%% CALCULATE X SPEs
x_hat = xs*xl';
error = x_fit_centred -x_hat;

SPE = sum(error.^2,2);

m = mean(SPE);
v = var(SPE);

h = 2*m^2/v;
g = v/(2*m);
chi_95 = g*chi2inv(.95, h);
chi_99 = g*chi2inv(.99, h);

%%%% Save Resulting Model %%%%



Q_model_Final.coeffs = beta_stored(:,:,nPC_idx);
Q_model_Final.y_index = desired_outputs;
Q_model_Final.x_index = [];
Q_model_Final.x0_index = desired_initial_states;
Q_model_Final.u_index = desired_inputs;
Q_model_Final.RMSE = RMSE;

Q_model_Final.f_index = [];
Q_model_Final.q_index = 1;
Q_model_Final.q_times = [];
Q_model_Final.nb_fit = nb_fit;
Q_model_Final.nb_val = nb_val;
Q_model_Final.nb_alt = (nb_test+nb_fit+nb_val) - nb_fit - nb_val;
Q_model_Final.nPC = PCs_to_try(nPC_idx);
Q_model_Final.sa = std(xs);
Q_model_Final.T95 = (nb_fit-1)*(nb_fit+1)*Q_model_Final.nPC/(nb_fit*(nb_fit-Q_model_Final.nPC))*finv(.95,Q_model_Final.nPC,nb_fit-Q_model_Final.nPC);
Q_model_Final.T99 = (nb_fit-1)*(nb_fit+1)*Q_model_Final.nPC/(nb_fit*(nb_fit-Q_model_Final.nPC))*finv(.99,Q_model_Final.nPC,nb_fit-Q_model_Final.nPC);
Q_model_Final.wstar = xl;
Q_model_Final.wstarAll = wstarAll;
Q_model_Final.Q = yl;
Q_model_Final.P = P;
Q_model_Final.W = W;
Q_model_Final.PAll = PAll;
Q_model_Final.OMEGA = Q_model_Final.wstar'*(x_fit_centred'*x_fit_centred)*Q_model_Final.wstar./(size(x_fit_centred,1)-1);

Q_model_Final.SPE95 = chi_95;
Q_model_Final.SPE99 = chi_99;

Q_model_Final.SPE_lim_calc = @(Prob) g*chi2inv(Prob, h);
Q_model_Final.T_lim_calc = @(Prob) (nb_fit-1)*(nb_fit+1)*Q_model_Final.nPC/(nb_fit*(nb_fit-Q_model_Final.nPC))*finv(Prob,Q_model_Final.nPC,nb_fit-Q_model_Final.nPC);

Q_model_Final.x_mean = x_fit_mu;
switch model_use
    case 'MultiMod'
        Q_model_Final.x_std = x_fit_sigma.*...
            [ones(1,size(Q_model_Final.x0_index,2))...
            ones(1,size(x_fit_sigma,2)-size(Q_model_Final.x0_index,2))*sqrt(size(x_fit_centred,2))];
    case 'MissingData'
        Q_model_Final.x_std = x_fit_sigma;
end
Q_model_Final.y_mean = qual_fit_mu;
Q_model_Final.y_std = qual_fit_sigma;



Q_model_Final.error = [0 0 0];

save(fname,'Q_model_Final');
%%
figure;

subplot(3,1,1);
plot(qual_test(:,1),'.');
hold on;
plot(qual_pred(:,1),'or');
hold off;
ylabel('Glucose g/L');
xlim([0,11]);

subplot(3,1,2);
plot(qual_test(:,2),'.');
hold on;
plot(qual_pred(:,2),'or');
hold off;
ylabel('Biomass g/L');
xlim([0,11]);

subplot(3,1,3);
plot(qual_test(:,3),'.');
hold on;
plot(qual_pred(:,3),'or');
hold off;
ylabel('Penicillin g/L');
xlim([0,11]);

legend('Actual','Predicted');
xlabel('Batch');

figure
scatter3(qual_test(:,1),qual_test(:,2),qual_test(:,3),'.');
hold on
scatter3(qual_pred(:,1),qual_pred(:,2),qual_pred(:,3),'or')
xlabel('Glucose g/L');
ylabel('Biomass g/L');
zlabel('Penicillin g/L');
legend('Actual','Predicted');