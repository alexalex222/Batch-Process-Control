clear;
close all;
clc;

% X_std = [.01 .07 .005 .001 4 .05 .005];
% Noise_std = [.01 .01 .001 .003 .01 .01]*0.01;
% bdb = BatchDatabase(size([15 1.16 0.1 0 100 0.5 298],2), 6, 6, 3, Tf, samt);
% bdb.AddICStd(X_std');
% bdb.AddMeasStd(Noise_std');
% ICs=abs(bdb.GetICs(randi(100), 7));
% test the Penicillin Simulator

t=0.5;
x_now= [15 1.16 0.1 0 100 0.5 298];
u_mpc=[0.045; 8; 30];
Terror_past=0;
jacket_temp_ideal_past=298;
[y_next, x_next, q_next, Terror_current, jacket_temp_ideal_current]=...
    PenicillinSimulator(t,x_now, u_mpc,Terror_past,jacket_temp_ideal_past);
X=zeros(601,7);
Y=zeros(601,6);
jacket_temp_ideal=zeros(601,1);
X(1,:)=x_now;
Y(1,:)=[x_now(7) x_now(7) x_now(2) x_now(6) x_now(5) 296];

for i= 2: 601
    x_now=x_next;
    X(i,:)=x_next;
    Y(i,:)=y_next;
    jacket_temp_ideal(i,:)=jacket_temp_ideal_current;
    Terror_past=Terror_current;
    jacket_temp_ideal_past=jacket_temp_ideal_current;
    [y_next, x_next,q_next, Terror_current, jacket_temp_ideal_current,noise]=...
    PenicillinSimulator(t,x_now, u_mpc,Terror_past,jacket_temp_ideal_past);
end

figure;
subplot 321
plot(Y(:,1));
title('Y1');

subplot 322
plot(Y(:,2));
title('Y2');

subplot 323
plot(Y(:,3));
title('Y3');

subplot 324
plot(Y(:,4));
title('Y4');

subplot 325
plot(Y(:,5));
title('Y5');

subplot 326
plot(Y(:,6));
title('Y6');


load test8

figure;
subplot 321
plot(bdb.data.ynom(1,:));
title('Y1');

subplot 322
plot(bdb.data.ynom(2,:));
title('Y2');

subplot 323
plot(bdb.data.ynom(3,:));
title('Y3');

subplot 324
plot(bdb.data.ynom(4,:));
title('Y4');

subplot 325
plot(bdb.data.ynom(5,:));
title('Y5');

subplot 326
plot(bdb.data.ynom(6,:));
title('Y6');
