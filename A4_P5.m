% Derek Woodard
% 100827083
% Assignment 4

close all
clear all
set(0,'DefaultFigureWindowStyle','docked');

%% Part 5
R1 = 1;
g1 = 1/R1;
C1 = 0.25;
R2 = 2;
g2 = 1/R2;
L = 0.2;
alpha = 100;
R4 = 0.1;
g4 = 1/R4;
Ro = 1000;
go = 1/Ro;
Cn = 0.00001;
In = 0.001;

nodes = 5;

global G C b;

% calculate R3 by going back to A3 a sweeping voltage across the bottleneck
% from 0.1-10V to extract a resistance value
% The value is calculated in A4_P1_using_A3.m
% The Resistance value calculated is always around 184
% Therefore, R3 = 184 Ohms is used
R3 = 184;
% R3 = 15;

%%
count = 1;
time = 1;
steps = 1000;
dt = time/steps;
times = zeros(steps,1);
Vouts_step = zeros(steps,1);
Vins_step = zeros(steps,1);

%% ci

if i*dt < 0.03
    vin = 0;
else
    vin = 1;
end

% Using the above KCL equations, the G and C matrix are formed
G = sparse(nodes, nodes);
C = sparse(nodes, nodes);
b = sparse(nodes,1);

% Using the stamps written in ELEC 4506, the G and C matricies are
% formed
res(1,2,R1);
res(2,0,R2);
res(3,0,R3);
res(4,5,R4);
res(5,0,Ro);
cap(1,2,C1);
cap(3,0,Cn);
ind(2,3,L);
vcvs(4,0,3,0,alpha/R3);
vol(1,0,vin);
cur(3,0,In);

if i == 1
    A = G + C./dt;
    tempb = b;
else
    A = G + C./dt;
    tempb = b + (C./dt)*prevV;
end

[Low,Up,P,Q] = lu(A,0.1);
z = Low\(P*tempb);
y = Up\z;
V = Q*y;   
prevV = V;

Vins_step(count) = vin;
Vouts_step(count) = V(5);
times(count) = i*dt;
count = count+1;


