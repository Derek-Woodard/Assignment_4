% Derek Woodard
% 100827083
% Assignment 4

close all
clear all
set(0,'DefaultFigureWindowStyle','docked');

%% Part 1
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
% Following KCL, the following equations are found
% (1) V1 = Vin
% (2) G1*(V1-V2) + C*((d(V1-V2))/dt) + I1 = 0
% (3) G1*(V2-V1) + C*((d(V2-V1))/dt) + G2*V2 + IL = 0
% (4) V2 - V3 - L*(d(I1)/dt) = 0
% (5) G3*V3 - IL = 0
% (6) G4*(V5 - V4) + Go*V5
% (7) G4*(V4 - V5) + Iv = 0
% (8) V4 - alpha*G3*V3 = 0

%%
% The KCL equations in the frequency domain
% (1) V1 = Vin
% (2) G1*(V1-V2) + C*(j*omega*(V1-V2)) + I1 = 0
% (3) G1*(V2-V1) + C*(j*omega*(V2-V1)) + G2*V2 + IL = 0
% (4) V2 - V3 - L*(j*omega*I1) = 0
% (5) G3*V3 - IL = 0
% (6) G4*(V5 - V4) + Go*V5
% (7) G4*(V4 - V5) + Iv = 0
% (8) V4 - alpha*G3*V3 = 0

Vins = zeros(80,1);
VoutsDC = zeros(80,1);
V3s = zeros(80,1);
count = 1;

%% 3i
% Now we sweep the input voltage from -10 to 10V
for vin = -10:0.25:10
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
    ind(2,3,L);
    vcvs(4,0,3,0,alpha/R3);
    vol(1,0,vin);
    
    if count == 1
        [Low,Up,P,Q] = lu(G,0.1);
    end
    z = Low\(P*b);
    y = Up\z;
    V = Q*y;
    
    Vins(count) = vin;
    VoutsDC(count) = V(5);
    V3s(count) = V(3);
    count = count+1;
end

figure(1)
plot(Vins, VoutsDC)
hold on
plot(Vins, V3s)
title('DC voltage sweep')
xlabel('Vin (V)');
ylabel('Voltage (V)')
legend('Vout', 'V3')

%% 3ii
% Now we plot the AC case
VoutsAC = zeros(500,1);
omegas = zeros(500,1);
count = 1;

for freq = 0:0.1:100
    omega = 2*pi*freq;
    omegas(count) = omega;
    A = G + j*omega*C;

    [Low,Up,P,Q] = lu(A,0.1);
    z = Low\(P*b);
    y = Up\z;
    V = Q*y;
    
    VoutsAC(count) = V(5);
    count = count+1;
end

figure(2)
plot(omegas, VoutsAC)
title('Vout as a function of omega')
xlabel('omega (rads/s)')
ylabel('Vout (V)')

gains = 20.*log(VoutsAC./vin);
figure(3)
plot(omegas,gains)
title('Gain vs. omega')
xlabel('omega (rads/s)')
ylabel('Gain (dB)')

%% 3iii
% Now we plot gain with random perturbations on C
omega = pi;
std = 0.5;
gainspert = zeros(500,1);

for i = 1:499
    Cpert = C;
    Cpert(1:2,1:2) = Cpert(1:2,1:2)+(0.05*randn);
    A = G + j*omega*Cpert;
    
    [Low,Up,P,Q] = lu(A,0.1);
    z = Low\(P*b);
    y = Up\z;
    V = Q*y;
    
    gainspert(i) = 20.*log(abs(V(5))/vin);
end

figure(4)
histogram(gainspert)
title('Gain with perturbations on C around 0.05')
xlabel('Gain (dB)')
ylabel('Instances')