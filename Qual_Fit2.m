clear variables;
close all;


load 'test4' 

fname = ['Q_model_V1'];

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


%get data from database and initial conditions

ynom = bdb.data.ynom(:,:);
unom = bdb.data.unom(:,:);
initialcond_nom = bdb.data.xnom(:,1);

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

% I believe it is here that the code works for you to say the outputs that you want to put in
% the model
%concatenates by adding the desired outputs, with the final row as input
% index = 1;
%for i=[1,4,5,6]
% for i=[1,2,5,6]
%     if i<6
%         concat_fit(index,:,:) = y_fit(i,:,:);
%         concat_val(index,:,:) = y_val(i,:,:);
%         concat_test(index,:,:) = y_test(i,:,:);
%         concat_nom(index,:) = ynom(i,:);
%     else
%         if(use_input)
%             concat_fit(index,:,:) = u_fit(1,:,:);
%             concat_val(index,:,:) = u_val(1,:,:);
%             concat_test(index,:,:) = u_test(1,:,:);
%             concat_nom(index,:) = unom(1,:);
%         end
%     end
%     index = index + 1;
% end

% Why add only the first input ? Why add the input in the final row ?
%concatenates by adding the desired outputs, with the final row as input
% index = 1;
% for i=1:numel(desired_outputs)
%     if i<numel(desired_outputs)
%         concat_fit(index,:,:) = y_fit(desired_outputs(i),:,:);
%         concat_val(index,:,:) = y_val(desired_outputs(i),:,:);
%         concat_test(index,:,:) = y_test(desired_outputs(i),:,:);
%         concat_nom(index,:) = ynom(desired_outputs(i),:);
%     else
%         if(use_input)
%             concat_fit(index,:,:) = u_fit(1,:,:);
%             concat_val(index,:,:) = u_val(1,:,:);
%             concat_test(index,:,:) = u_test(1,:,:);
%             concat_nom(index,:) = unom(1,:);
%         end
%     end
%     index = index + 1;
% end


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
yvals = sqrt((T_sqrd_99_2comp-xvals.^2./s1^2)*s2.^2);

hold on
plot(xvals, yvals, '--k');
yvals = -yvals;
plot(xvals, yvals, '--k');

% subplot(1,2,2);
% scatter(ys(:,1),ys(:,2))

% figure
% subplot(2,2,1);
% bar(xl(1:3,2))
% subplot(2,2,2);
% bar(xl(1:3,1))
% subplot(2,2,3);
% bar(yl(:,1))

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
% 
% 
% 
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



Q_model.coeffs = beta_stored(:,:,nPC_idx);
Q_model.y_index = desired_outputs;
Q_model.x_index = [];
Q_model.x0_index = desired_initial_states;
Q_model.u_index = desired_inputs;
Q_model.RMSE = RMSE;

Q_model.f_index = [];
Q_model.q_index = 1;
Q_model.q_times = [];
Q_model.nb_fit = nb_fit;
Q_model.nb_val = nb_val;
Q_model.nb_alt = (nb_test+nb_fit+nb_val) - nb_fit - nb_val;
Q_model.nPC = PCs_to_try(nPC_idx);
Q_model.sa = std(xs);
Q_model.T95 = (nb_fit-1)*(nb_fit+1)*Q_model.nPC/(nb_fit*(nb_fit-Q_model.nPC))*finv(.95,Q_model.nPC,nb_fit-Q_model.nPC);
Q_model.T99 = (nb_fit-1)*(nb_fit+1)*Q_model.nPC/(nb_fit*(nb_fit-Q_model.nPC))*finv(.99,Q_model.nPC,nb_fit-Q_model.nPC);
Q_model.wstar = xl;
Q_model.wstarAll = wstarAll;
Q_model.Q = yl;
Q_model.P = P;
Q_model.W = W;
Q_model.PAll = PAll;
Q_model.OMEGA = Q_model.wstar'*(x_fit_centred'*x_fit_centred)*Q_model.wstar./(size(x_fit_centred,1)-1);

Q_model.SPE95 = chi_95;
Q_model.SPE99 = chi_99;

Q_model.SPE_lim_calc = @(Prob) g*chi2inv(Prob, h);
Q_model.T_lim_calc = @(Prob) (nb_fit-1)*(nb_fit+1)*Q_model.nPC/(nb_fit*(nb_fit-Q_model.nPC))*finv(Prob,Q_model.nPC,nb_fit-Q_model.nPC);

Q_model.x_mean = x_fit_mu;
switch model_use
    case 'MultiMod'
        Q_model.x_std = x_fit_sigma.*...
            [ones(1,size(Q_model.x0_index,2))...
            ones(1,size(x_fit_sigma,2)-size(Q_model.x0_index,2))*sqrt(size(x_fit_centred,2))];
    case 'MissingData'
        Q_model.x_std = x_fit_sigma;
end
Q_model.y_mean = qual_fit_mu;
Q_model.y_std = qual_fit_sigma;



Q_model.error = [0 0 0];

save(fname,'Q_model');

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

% figure
% 
% subplot(2,2,1)
% scatter(qual_test(:,1),qual_pred(:,1))
% % axis([0.7 0.76 0.7 0.76])
% % line([0.7 0.76],[0.7 0.76])
% xlabel('Actual Glucose')
% ylabel('Predicted Glucose')
% 
% subplot(2,2,2)
% scatter(qual_test(:,2),qual_pred(:,2))
% % axis([86000 106000 86000 106000])
% % line([86000 106000],[86000 106000])
% xlabel('Actual Biomass')
% ylabel('Predicted Biomass')
% 
% subplot(2,2,3)
% scatter(qual_test(:,3),qual_pred(:,3))
% % axis([150000 185000 150000 185000])
% % line([150000 185000],[150000 185000])
% xlabel('Actual Pen')
% ylabel('Predicted Pen')

% qual_pred_centred = x_test_centred*beta_stored(:,:,3);
% for i=1:size(qual_pred_centred,1)
%     for j = 1:size(qual_pred_centred,2)
%         qual_pred(i,j) = qual_pred_centred(i,j)*qual_fit_sigma(j)+qual_fit_mu(j);
%         error(i,j) = abs(qual_pred(i,j) - qual_test(i,j));
%     end
% end
% 
% 
% subplot(3,3,4)
% scatter(qual_test(:,1),qual_pred(:,1))
% % axis([0.7 0.76 0.7 0.76])
% % line([0.7 0.76],[0.7 0.76])
% xlabel('Actual Conversion')
% ylabel('Predicted Conversion')
% 
% subplot(3,3,5)
% scatter(qual_test(:,2),qual_pred(:,2))
% % axis([86000 106000 86000 106000])
% % line([86000 106000],[86000 106000])
% xlabel('Actual NAvg MW')
% ylabel('Predicted NAvg MW')
% 
% subplot(3,3,6)
% scatter(qual_test(:,3),qual_pred(:,3))
% % axis([150000 185000 150000 185000])
% % line([150000 185000],[150000 185000])
% xlabel('Actual WAvg MWw')
% ylabel('Predicted WAvg MWw')
% 
% qual_pred_centred = x_test_centred*beta_stored(:,:,nPC_idx);
% for i=1:size(qual_pred_centred,1)
%     for j = 1:size(qual_pred_centred,2)
%         qual_pred(i,j) = qual_pred_centred(i,j)*qual_fit_sigma(j)+qual_fit_mu(j);
%         error(i,j) = abs(qual_pred(i,j) - qual_test(i,j));
%     end
% end
% 
% subplot(3,3,7)
% scatter(qual_test(:,1),qual_pred(:,1))
% % axis([0.7 0.76 0.7 0.76])
% % line([0.7 0.76],[0.7 0.76])
% xlabel('Actual Conversion')
% ylabel('Predicted Conversion')
% 
% subplot(3,3,8)
% scatter(qual_test(:,2),qual_pred(:,2))
% % axis([86000 106000 86000 106000])
% % line([86000 106000],[86000 106000])
% xlabel('Actual NAvg MW')
% ylabel('Predicted NAvg MW')
% 
% subplot(3,3,9)
% scatter(qual_test(:,3),qual_pred(:,3))
% % axis([150000 185000 150000 185000])
% % line([150000 185000],[150000 185000])
% xlabel('Actual WAvg MWw')
% ylabel('Predicted WAvg MWw')

% figure
% hist(error(:,1),50)