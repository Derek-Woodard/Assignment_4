% Derek Woodard
% 100827083
% Assignment 4

close all
clear all
set(0,'DefaultFigureWindowStyle','docked');

%% Part 4
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
count = 1;
time = 1;
steps = 1000;
% steps = 2000;
% steps = 5000;
dt = time/steps;
times = zeros(steps,1);
Vouts_step = zeros(steps,1);
Vins_step = zeros(steps,1);

%% di
for i = 1:steps
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
    ind(2,3,L);
    vcvs(4,0,3,0,alpha/R3);
    vol(1,0,vin);

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
end

figure(1)
subplot(2,1,1)
plot(times,Vins_step)
hold on
plot(times,Vouts_step)
title('In & Out voltage with step input')
xlabel('time (s)')
ylabel('Voltage (V)')
legend('Vin', 'Vout')

%% dii
freq = 1/0.03;
count = 1;
times = zeros(steps,1);
Vouts_sin = zeros(steps,1);
Vins_sin = zeros(steps,1);
for i = 1:steps
    vin = sin(2*pi*freq*(i*dt));
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
    
    Vins_sin(count) = vin;
    Vouts_sin(count) = V(5);
    times(count) = i*dt;
    count = count+1;
end

figure(2)
subplot(2,1,1)
plot(times,Vins_sin)
hold on
plot(times,Vouts_sin)
title('In & Out voltage with sin input')
xlabel('time (s)')
ylabel('Voltage (V)')
legend('Vin', 'Vout')

%% diii
count = 1;
times = zeros(steps,1);
Vouts_gaus = zeros(steps,1);
Vins_gaus = zeros(steps,1);

prange = 1:1:180;
gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
pulse = gaus(prange,90,30,1,0);
Vins_gaus(1:150) = pulse(31:180)';

for i = 1:steps
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
    vol(1,0,Vins_gaus(count));

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

figure(3)
subplot(2,1,1)
plot(times,Vins_gaus)
hold on
plot(times,Vouts_gaus)
title('In & Out voltage with guassian pulse input')
xlabel('time (s)')
ylabel('Voltage (V)')
legend('Vin', 'Vout')

%%
% Now we plot the frequency content of the input and output signals
freqrange_step = (-length(Vins_step)/2:length(Vins_step)/2-1);
figure(1)
subplot(2,1,2)
hold on
plot(freqrange_step,mag2db(abs(fftshift(fft(Vins_step)))));
plot(freqrange_step,mag2db(abs(fftshift(fft(Vouts_step)))));
hold off
title('Frequency Response')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
legend('Vin', 'Vout')

freqrange_sin = (-length(Vins_sin)/2:length(Vins_sin)/2-1);
figure(2)
subplot(2,1,2)
hold on
plot(freqrange_sin,mag2db(abs(fftshift(fft(Vins_sin)))));
plot(freqrange_sin,mag2db(abs(fftshift(fft(Vouts_sin)))));
hold off
title('Frequency Response')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
legend('Vin', 'Vout')

freqrange_gaus = (-length(Vins_gaus)/2:length(Vins_gaus)/2-1);
figure(3)
subplot(2,1,2)
hold on
plot(freqrange_gaus,mag2db(abs(fftshift(fft(Vins_gaus)))));
plot(freqrange_gaus,mag2db(abs(fftshift(fft(Vouts_gaus)))));
hold off
title('Frequency Response')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
legend('Vin', 'Vout')

