close all
clear all
clc

%Some parameters in ODEs
a=1.16; %a is Dissolved Oxygen Concentrarion at saturation in g/L. In the article, this variable is cL and it is explict in equations (2), (3), (10) and (12). This value comes from pensim line 23.
F=0; %F is Feed Flow Rate of Substrate in L/h. In the article, this variable is also named F and it is explicit in equations (5), (11), (14) and (17). This value comes from pensim2 line 142.
fgg=8; %fgg is Flow Rate of Oxygen. In the article, this variable is named fg and it is explicit in equation (13). This value comes from pensim line 34.
pww=30; %pww is Power Input. In the article, this variable is named Pw and it is explicit in equation (13). This value comes from pensim line 34.
Tf=296; % Tf is Feed Temperature of Substrate in K. In the article, this variable is named Tf and it is explicit in equation (17). This value comes from pensim line 34.
Tset = 298; % Set point for the batch temperature [K]
Fc=.06; %Fc is Cooling Water Flow Rate in L/h. In the article, this variable is also named Fc and it is explicit in equation (17). This value comes from pensim line 149.
b=0; %b => flag for heating/cooling switching (b=0 : Cooling, b=1 : Heating). It is explicit in equation (17). This value comes from pensim2 line 159.
tfl=1; %tfl is a variable used for controller that controls the temperature. If tfl=0 means open loop for temperature. This value comes from pensim2 line 259.
pH=5; %This is pH that will be used to calculate the concentration of the hydrogen ion. This value comes from pH Control section ( 2.3 ) in the paper.
ti=0; %Initial time of integration.
t_feed_start = 40; % The time when a feed will be started
F_feed = 4; %0.0450; % Feed rate for fed batch part
tf=300; %Final time of integration. This value comes from pensim line 72.
samt = .5; % sampling time 

nFE = tf/samt; % number of sampling instances

 inp1=[a F fgg pww Tf Fc 300 pH]; %Inputs to the ODE file pengsimv2_0
 % Proposed input <
 %inp1=[a F fgg pww Tf Fc b tfl pH F_feed]; %Inputs to the ODE file pengsimv2_0

 Y0 = [15 1.16 0.1 0 100 0.5 298]; %Initial conditions for the states.
% Y0(1) Initial substrate concentration
% Y0(2) Initial Dissolved O2
% Y0(3) Initial Biomass concentration
% Y0(4) Initial penicillin concentration
% Y0(5) Initial fermentor volume
% Y0(6) Initial CO2
% Y0(7) Initial fermentor temperature
opts = odeset('reltol', 1e-9);

X = zeros(nFE+1,size(Y0,2)); % allocate memory for x
Terror = zeros(nFE,1);
u = zeros(nFE, 4); % Glucose Feed, Aeration, Agitator power, Temp of coil
 u_implemented = zeros(nFE,1);

% Best so far:
% Kc = -100;
% tauI = 20000;

Kc = -4000;
tauI = 1000;

% Kc =-1;
% tauI = 1;

X(1,:) = Y0; % initialize the states


%Solve ODEs using ode23s
for k = 1:nFE
    Terror(k) = (X(k,7)-Tset);
    if(k == 1)
        u(k, 4) = Kc*Terror(k) + 298;
    else
        delta_u = Kc*(Terror(k) - Terror(k-1) + samt/tauI*Terror(k));
        u(k, 4) = u(k-1,4) + delta_u;
    end
    u_implemented(k) = u(k,4);
    
    if(u(k,4)<290)
       u_implemented(k) = 290; 
    elseif(u(k,4)>323)
        u_implemented(k) = 323;
    end
    
    inp1(7) = u_implemented(k);
    
    [T_out,Y_out]=ode15s('pengsimv2_1',[(k-1) k]*samt,X(k,:),opts,inp1);
    
    X(k+1,:) = Y_out(end,:);

end

fprintf(['SSE in temp control was: ' num2str(sum(Terror(:).^2)) '\n']);
figure
%Substrate Concentration
subplot(3,3,1)
plot(T1,Y1(:,1)) %Plot Substrate Concentration
xlabel('Time (h)') %x axis
ylabel('Substrate Concentration (g/L)') %y axis
%Concentration of Oxygen Dissolved
subplot(3,3,2)
plot(T1,Y1(:,2)) %Plot Concentration of Oxygen Dissolved
xlabel('Time (h)') %x axis
ylabel('Oxygen Dissolved (g/L)') %y axis
%Biomass Concentration
subplot(3,3,3)
plot(T1,Y1(:,3)) %Plot Biomass Concentration
xlabel('Time (h)') %x axis
ylabel('Biomass Concentration (g/L)') %y axis
%Penicillin Concentration
subplot(3,3,4)
plot(T1,Y1(:,4)) %Plot Penicillin Concentration
xlabel('Time (h)') %x axis
ylabel('Penicillin Concentration (g/L)') %y axis
%Volume
subplot(3,3,5)
plot(T1,Y1(:,5)) %Volume
xlabel('Time (h)') %x axis
ylabel('Fermentor Volume (L)') %y axis
%CO2 Evolution
subplot(3,3,6)
plot(T1,Y1(:,6)) %CO2 Evolution
xlabel('Time (h)') %x axis
ylabel('Concentration of CO2 (mmolCO2/L)') %y axis
%Temperature
subplot(3,3,7)
plot(T1,Y1(:,7)) %Temperature
xlabel('Time (h)') %x axis
ylabel('Temperature (K)') %y axis



%Solve ODEs using ode15s

[T2,Y2]=ode15s('pengsimv2_0',[ti t_feed_start],Y0,opts,inp1);
inp1(2) = F_feed;

[T3,Y3]=ode15s('pengsimv2_0',[t_feed_start tf],Y2(end,:),opts,inp1);

T2 = [T2; T3];
Y2 = [Y2; Y3];

figure
%Substrate Concentration
subplot(3,3,1)
plot(T2,Y2(:,1)) %Plot Substrate Concentration
xlabel('Time (h)') %x axis
ylabel('Substrate Concentration (g/L)') %y axis
%Concentration of Oxygen Dissolved
subplot(3,3,2)
plot(T2,Y2(:,2)) %Plot Concentration of Oxygen Dissolved
xlabel('Time (h)') %x axis
ylabel('Oxygen Dissolved (g/L)') %y axis
%Biomass Concentration
subplot(3,3,3)
plot(T2,Y2(:,3)) %Plot Biomass Concentration
xlabel('Time (h)') %x axis
ylabel('Biomass Concentration (g/L)') %y axis
%Penicillin Concentration
subplot(3,3,4)
plot(T2,Y2(:,4)) %Plot Penicillin Concentration
xlabel('Time (h)') %x axis
ylabel('Penicillin Concentration (g/L)') %y axis
%Volume
subplot(3,3,5)
plot(T2,Y2(:,5)) %Volume
xlabel('Time (h)') %x axis
ylabel('Fermentor Volume (L)') %y axis
%CO2 Evolution
subplot(3,3,6)
plot(T2,Y2(:,6)) %CO2 Evolution
xlabel('Time (h)') %x axis
ylabel('Concentration of CO2 (mmolCO2/L)') %y axis
%Temperature
subplot(3,3,7)
plot(T2,Y2(:,7)) %Temperature
xlabel('Time (h)') %x axis
ylabel('Temperature (K)') %y axis


