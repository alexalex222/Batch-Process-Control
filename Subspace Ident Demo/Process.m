clear variables;
close all;

k_1 = .03;
k_2 = .02;

CA_in = 2.0;
CB_in = .5;
V = 20;
v = 4;

% A -> B -> C

u0 = [CA_in; CB_in];

A_cont = [ -k_1-v/V 0 0;
    k_1 -k_2-v/V 0;
    0 k_2 -v/V];

B_cont = [v/V 0; 0 v/V; 0 0];

C_cont = [ 0 1 1];

xss = A_cont\(-B_cont*u0);