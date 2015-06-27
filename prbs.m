% Input PRBS Sequence design
%
% I/O structure : [u,z]=prbs(sp,tt,step2,aa)
%
% Outputs : u => PRBS,  z => Filtered PRBS
% Inputs : sp => Set points (row vector), tt : Process time span, 
%          step2 => Sampling time, aa => Set point variations (row vector)

% Example => sp=[8 30 0.045 296 5 298];
%				aa=[0.3 1 0.001 0.5 0.03 0.03];
%
% Ref. : System Identification, Theory for the user
% by L. Ljung, Prentice Hall, 1987, page 373.
%
% written by Cenk Undey (c) 01/26/2000, IIT, CHEE
%
% All rights reserved. Copyright (c) IIT 2002.

function z=prbs(sp,tt,step2,aa)

u=[];z=[];
if nargin==3
    aa=[0.3 1 0.005 0.5 0.03 0.03];% Deviations from set points
end

u1=sp+aa;
u2=sp-aa;

fr=ceil(20/step2); % Frequency

b=ones(1,fr)/fr;

r=sp/0.01;
t=ceil(tt/step2);

for g=1:length(sp)
   uu=sp(g);
   w=randn(1,t);
   for i=1:t
      uu(i+1)=0.5*(u1(g)+u2(g))+0.5*(u1(g)-u2(g))*sign(uu(i)/r(g)-w(i));
   end
   
   uu(length(uu))=sp(g);%PRBS
   
   zz=filtfilt(b,1,uu);%Filtered PRBS Matlba's Signal Processing Toolbox function filtfilt.m
   
   z=[z; zz];u=[u; uu];
end


   
   