% Derek Woodard
% 100827083
% Assignment 4

close all
clear all
set(0,'DefaultFigureWindowStyle','docked');

%% Part 6
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
% A capacitor and current source are added in parallel with R3 to simulate
% noise
time = 1;
steps = 1000;
dt = time/steps;

prange = 1:1:180;
gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
pulse = gaus(prange,90,30,1,0);

Vc = 1;
Vins_gaus = zeros(steps,1);
Vins_gaus(1:150) = pulse(31:180)';

count = 1;

times = zeros(steps,1);

Vouts_gaus = zeros(steps,1);

gaus_In = In * randn(1, length(Vins_gaus));
for i = 1:steps
    % Using the above KCL equations, the G, C, and b matricies are formed
    G = sparse(nodes, nodes);
    C = sparse(nodes, nodes);
    b = sparse(nodes,1);

    % Using the stamps written in ELEC 4506, the G, C, and b matricies are
    % formed
    res(1,2,R1);
    res(2,0,R2);
    res(3,0,R3);
    res(4,5,R4);
    res(5,0,Ro);
    cap(1,2,C1);
    cap(3,0,Cn(Vc));
    ind(2,3,L);
    vcvs(4,0,3,0,alpha/R3);
    vol(1,0,Vins_gaus(i));
    cur(3,0,gaus_In(i));

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

    Vouts_gaus(count) = V(5);
    times(count) = i*dt;
    count = count+1;
end

figure(1)
subplot(2,1,1)
plot(times,Vins_gaus)
hold on
plot(times,Vouts_gaus)
title('In & Out Voltage with Guassian Pulse Input with Noise')
xlabel('Time (s)')
ylabel('Voltage (V)')
legend('Vin', 'Vout')

freqrange_gaus = (-length(Vins_gaus)/2:length(Vins_gaus)/2-1);
subplot(2,1,2)
hold on
plot(freqrange_gaus,mag2db(abs(fftshift(fft(Vins_gaus)))));
plot(freqrange_gaus,mag2db(abs(fftshift(fft(Vouts_gaus)))));
hold off
title('Frequency Response')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
legend('Vin', 'Vout')


figure(2)
histogram(gaus_In);
title('Gaussian Distribution of Random In Values')
xlabel('In (A)')
ylabel('Instances')

