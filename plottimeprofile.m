clc;
clear;
close all;

load test4;

X = bdb.data.xnom;
X = X';
X = zscore(X);

bias = repmat(min(X),601,1);
X = X-bias;
X(:,2) = X(:,2)./2.7+0.5;
X(:,4) = X(:,4).*1.5;

figure
plot(X(:,1),'LineWidth',3);
hold on;
plot(X(:,2),'r','LineWidth',3);
plot(X(:,3),'g','LineWidth',3);
plot(X(:,4),'k','LineWidth',3);
xlim([0 600]);
legend('Glucose','Dissolved O_2','Biomass','Penicillin');
xlabel('Time');
ylabel('Concentration');
set(gca,'XTick',[]);
set(gca,'YTick',[]);

figure;
subplot 221
plot(bdb.data.y(2,:,1));
title('Reactor Temperature');
xlim([0 600]);


subplot 222
plot(bdb.data.y(3,:,1));
title('Dissolved O2');
xlim([0 600]);


subplot 223
plot(bdb.data.y(4,:,1));
title('Dissoved CO2');
xlim([0 600]);


subplot 224
plot(bdb.data.y(5,:,3));
title('Fermentor Volume');
xlim([0 600]);
