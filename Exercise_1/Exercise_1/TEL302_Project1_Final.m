%%%%%%            TEL302 - Project 1                %%%%%%
%%%%%% Ilektra-Despoina Papamatthaiaki (2018030106) %%%%%%

clear all;
close all;
clc;

% A1 

T = 10^(-3);
over = 10;
Ts = T / over;
A = 4;
a1 = 0;
a2 = 0.5;
a3 = 1;

[phi1, t1] = srrc_pulse(T, over, A, a1);
[phi2, t2] = srrc_pulse(T, over, A, a2);
[phi3, t3] = srrc_pulse(T, over, A, a3);

figure(1);
plot(t1, phi1, 'red');
hold on;
plot(t2, phi2, 'blue');
hold on;
plot(t3, phi3, 'green');
set(gcf,'color', 'w');

% A2 

Nf = 1024;
Fs = 1/Ts;
Fx = -Fs/2:Fs/Nf:Fs/2-Fs/Nf;

S1 = fftshift(fft(phi1,Nf)* Ts);
S2 = fftshift(fft(phi2,Nf)* Ts);
S3 = fftshift(fft(phi3,Nf)* Ts);

figure(2);
plot(Fx, abs(S1).^2 , 'red');
hold on;
plot(Fx, abs(S2).^2 , 'blue');
hold on;
plot(Fx, abs(S3).^2 , 'green');
set(gcf,'color', 'w');

figure(3);
semilogy(Fx, abs(S1).^2 , 'red');
hold on;
semilogy(Fx, abs(S2).^2 , 'blue');
hold on;
semilogy(Fx, abs(S3).^2 , 'green');
set(gcf,'color', 'w');

% A3 

BW1 = (1 + a1) / (2 * T)
BW2 = (1 + a2) / (2 * T)
BW3 = (1 + a3) / (2 * T)

c1 = T ./ 10^(3); 
c2 = T ./ 10^(5); 

figure(4);
semilogy(Fx, abs(S1).^2 , 'red');
hold on;
semilogy(Fx, abs(S2).^2 , 'blue');
hold on;
semilogy(Fx, abs(S3).^2 , 'green');
hold on;
y = c1 * (Fx > -Fs/2);
plot(Fx,y,'magenta');
hold on;
y1 = c2 * (Fx > -Fs/2);
plot(Fx,y1,'magenta');
set(gcf,'color', 'w');

% B1

i = 5;

for k = 0:3
    figure(i);
    
    
    t1new = t1(1):Ts:t1(end) + T * k;
    t2new = t2(1):Ts:t2(end) + T * k;
    t3new = t3(1):Ts:t3(end) + T * k;
    
    oldPh1Shift1 = [phi1   zeros(1, k * (T/Ts))];  % old signal with padding
    oldPh1Shift2 = [phi2   zeros(1, k * (T/Ts))];  
    oldPh1Shift3 = [phi3   zeros(1, k * (T/Ts))];  
    
    newPh1Shift1 = [zeros(1, k * (T/Ts)) phi1];    % new signal with padding 
    newPh1Shift2 = [zeros(1, k * (T/Ts)) phi2];     
    newPh1Shift3 = [zeros(1, k * (T/Ts)) phi3];     
    
    multiPhi1 = oldPh1Shift1 .* newPh1Shift1;      % multiplication of the 2 signals
    multiPhi2 = oldPh1Shift2 .* newPh1Shift2;      
    multiPhi3 = oldPh1Shift3 .* newPh1Shift3;      
    
    subplot(3,2,1)
    plot(t1 + T * k, phi1, 'magenta');
    hold on;
    plot(t1, phi1, 'blue');
    title(sprintf('SRRC Pulses with a=0 for k=%d',k));
    
    subplot(3,2,3)
    plot(t2 + T * k, phi2, 'magenta');
    hold on;
    plot(t2, phi2, 'blue');
    title(sprintf('SRRC Pulses with a=0.5 for k=%d',k));
    
    subplot(3,2,5)
    plot(t3 + T * k, phi3, 'magenta');
    hold on;
    plot(t3, phi3, 'blue');
    title(sprintf('SRRC Pulses with a=1 for k=%d',k));
    
    subplot(3,2,2)
    plot(t1new, multiPhi1, 'magenta');
    title(sprintf('Mult of SRRC Pulses with a=0 for k=%d',k));
    
    subplot(3,2,4)
    plot(t2new, multiPhi2, 'magenta');
    title(sprintf('Mult of SRRC Pulses with a=0.5 for k=%d',k));
    
    subplot(3,2,6)
    plot(t3new, multiPhi3, 'magenta');
    title(sprintf('Mult of SRRC Pulses with a=1 for k=%d',k));
    set(gcf,'color', 'w');
    
    Sum1 = sum(multiPhi1 .* Ts)
    Sum2 = sum(multiPhi2 .* Ts)
    Sum3 = sum(multiPhi3 .* Ts)
    
    i = i+1;

end

% C1

T = 10^-3;
over = 10;
a = 0.5;
A = 4;
N = 100;
Ts = T / over;

b = (sign(randn(N, 1)) + 1)/2;

% C2

X = bits_to_2PAM(b);

X_delta = 1/Ts * upsample(X, over); 
t = 0:Ts:N*Ts*over-Ts;

figure();
stem(t,X_delta);
set(gcf,'color', 'w');

% C3

[phi_new, t_new] = srrc_pulse(T, over, A, a);
X_conv = conv(X_delta,phi_new)*Ts;
t_conv = t(1) + t_new(1):Ts:t(end) + t_new(end);
figure();
plot(t_conv,X_conv);
set(gcf,'color', 'w');

phi_new_rev = phi_new(end:-1:1); %reflected signal
t_new_rev =  -t_new(end:-1:1);   %reflected time
t_conv2 = t_conv(1) + t_new_rev(1):Ts:t_conv(end) + t_new_rev(end);

Z = conv(X_conv,phi_new_rev)*Ts;

figure();
plot(t_conv2,Z);
set(gcf,'color', 'w');

figure();
plot(t_conv2,Z);
hold on;
stem([0:N-1]*T,X);
set(gcf,'color', 'w');



