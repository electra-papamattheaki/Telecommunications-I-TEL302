clear all;
close all;
clc;

%A1
T = 10^(-2);
over = 10;
Ts = T / over;
A = 4;
a = 0.5;

[phi, t] = srrc_pulse(T, over, A, a);

figure(1);
plot(t, phi, 'blue');
title('SRRC Pulse');
set(gcf,'color', 'w');

Nf = 2048;
Fs = 1/Ts;
Fx = -Fs/2:Fs/Nf:Fs/2-Fs/Nf;
S = fftshift(fft(phi,Nf)* Ts);

figure(2);
semilogy(Fx, abs(S).^2 , 'blue');
title('Power Spectral Density');
set(gcf,'color', 'w');

%A2
N = 100;
b = (sign(randn(N, 1)) + 1)/2;

Xn = bits_to_2PAM(b);
X_delta = 1/Ts*upsample(Xn,over);
t1= 0:Ts:N*Ts*over-Ts;

X = conv(X_delta,phi)*Ts; 
tX = t1(1)+t(1):Ts:t1(end)+t(end);

figure(3);
plot(tX,X);
set(gcf,'color', 'w');

%A3
Ttotal = tX(end) - tX(1);
S2 = fftshift(fft(X,Nf)* Ts);
Px = ((abs(S2)) .^ 2) / Ttotal;

% periodogramma me plot
figure(4);
plot(Fx, Px , 'blue');
set(gcf,'color', 'w');

% periodogramma me semilogy
figure(5);
semilogy(Fx, Px , 'blue');
set(gcf,'color', 'w');

Px_sum = 0;

for K = 0:500
    
    b = (sign(randn(N, 1)) + 1)/2;
    Xn = bits_to_2PAM(b);
    X_delta = 1/Ts*upsample(Xn,over);
    t1= 0:Ts:N*Ts*over-Ts;
    X = conv(X_delta,phi)*Ts; 
    tX = t1(1)+t(1):Ts:t1(end)+t(end);
    Ttotal = tX(end) - tX(1);
    S2 = fftshift(fft(X,Nf)* Ts);
    Px = ((abs(S2)) .^ 2) / Ttotal;
    Px_sum = Px_sum + Px;
end

Px_final = Px_sum / K;
figure(6);
semilogy(Fx, Px_final , 'blue');
set(gcf,'color', 'w');

Sx = (var(Xn)/T) .* (abs(S).^2);

figure(7);
semilogy(Fx, Px_final , 'blue');
hold on
semilogy(Fx, Sx , 'red');
title('Theoretical & Experimental Power Spectral Density');
set(gcf,'color', 'w');
legend('experimental','theoretical');

%A4
b2 = (sign(randn(N/2, 1)) + 1)/2;
Xn2 = bits_to_4PAM(b2);
X_delta2 = 1/Ts*upsample(Xn2,over);
t2= 0:Ts:(N/2)*Ts*over-Ts;

X2 = conv(X_delta2,phi)*Ts; %convolution of ä(t ? kT) with phi gives ?(t ? nT)
tX2 = t2(1)+t(1):Ts:t2(end)+t(end);

figure(8);
plot(tX2,X2);
set(gcf,'color', 'w');

Ttotal2 = tX2(end) - tX2(1);
S4 = fftshift(fft(X2,Nf)* Ts);
Px2 = ((abs(S4)) .^ 2) / Ttotal2;

figure(9);
plot(Fx, Px2 , 'blue');
set(gcf,'color', 'w');

figure(10);
semilogy(Fx, Px2 , 'blue');
set(gcf,'color', 'w');


Px_sum2 = 0;

for K = 0:500
    
    b2 = (sign(randn(N/2, 1)) + 1)/2;
    Xn2 = bits_to_2PAM(b2);
    X_delta2 = 1/Ts*upsample(Xn2,over);
    t2= 0:Ts:(N/2)*Ts*over-Ts;
    X2 = conv(X_delta2,phi)*Ts; 
    tX2 = t2(1)+t(1):Ts:t2(end)+t(end);
    Ttotal2 = tX2(end) - tX2(1);
    S4 = fftshift(fft(X2,Nf)* Ts);
    Px2 = ((abs(S4)) .^ 2) / Ttotal2;
    Px_sum2 = Px_sum2 + Px2;
end

Px_final2 = Px_sum2 / K;
figure(11);
semilogy(Fx, Px_final2 , 'blue');
set(gcf,'color', 'w');

Sx2 = (var(Xn2)/T) .* (abs(S).^2);

figure(12);
semilogy(Fx, Px_final2 , 'blue');
hold on
semilogy(Fx, Sx2 , 'red');
set(gcf,'color', 'w');
title('Theoretical & Experimental Power Spectral Density');
legend('experimental','theoretical');

%A5
T5 = 2*T;
over2 = 2*over;
Ts = T5 / over2;
Nf = 4096;
Fs = 1/Ts;
Fx = -Fs/2:Fs/Nf:Fs/2-Fs/Nf;

[phi5, t5] = srrc_pulse(T5, over2, A, a);
S5 = fftshift(fft(phi5,Nf)* Ts);

b5 = (sign(randn(N, 1)) + 1)/2;
Xn5 = bits_to_2PAM(b5);
X_delta5 = 1/Ts*upsample(Xn5,over2);
t6= 0:Ts:N*Ts*over2-Ts;

X5 = conv(X_delta5,phi5)*Ts;
tX5 = t6(1)+t5(1):Ts:t6(end)+t5(end);

Ttotal5 = tX5(end) - tX5(1);
S6 = fftshift(fft(X5,Nf)* Ts);
Px5 = ((abs(S6)) .^ 2) / Ttotal5;

figure(13);
plot(Fx,Px5)
set(gcf,'color', 'w');

figure(14);
semilogy(Fx,Px5)
set(gcf,'color', 'w');

Px_sum5 = 0;

for K = 0:500
    b5 = (sign(randn(N, 1)) + 1)/2;
    Xn5 = bits_to_2PAM(b5);
    X_delta5 = 1/Ts*upsample(Xn5,over2);
    t6= 0:Ts:N*Ts*over2-Ts;
    X5 = conv(X_delta5,phi5)*Ts; %convolution of ä(t ? kT) with phi gives ?(t ? nT)
    tX5 = t6(1)+t5(1):Ts:t6(end)+t5(end);
    Ttotal5 = tX5(end) - tX5(1);
    S6 = fftshift(fft(X5,Nf)* Ts);
    Px5 = ((abs(S6)) .^ 2) / Ttotal5;
    Px_sum5 = Px_sum5 + Px5;
end

Px_final5 = Px_sum5 / K(end);
figure(15);
semilogy(Fx, Px_final5 , 'blue');
set(gcf,'color', 'w');

Sx5 = (var(Xn5)/T5) .* (abs(S5).^2);

figure(16);
semilogy(Fx, Px_final5 , 'blue');
hold on
semilogy(Fx, Sx5 , 'red');
set(gcf,'color', 'w');
legend('experimental','theoretical');
title('Theoretical & Experimental Power Spectral Density');
set(gcf,'color', 'w');

%B1

F0 = 4; 
t_b = 1:0.01:pi; 

% plot 5 instances of the given function
figure(17)
for i=1:5
   X  = normrnd(0,1); 
   phi_b = 2*pi*rand(1); 
   Y = X*cos(2*pi*F0*t_b+phi_b);
   plot(t_b,Y)
   hold on;
end

title('5 Instances of Y(t) = X*cos(2*pi*F0*t+phi)')
grid on;
set(gcf,'color', 'w');

%B2

% Ts = T/over;
% [phi, t] = srrc_pulse(T, over, A, a);
% Nf = 2048;
% Fs = 1/Ts;
% 
% X  = normrnd(0,1); 
% phi_b = 2*pi*rand(1);
% 
% Fx = -Fs/2:Fs/Nf:Fs/2-Fs/Nf;
% S = fftshift(fft(phi,Nf)*Ts);
% t1 = 0:Ts:N*Ts*over-Ts;
% tX = t1(1)+t(1):Ts:t1(end)+t(end);
% fo = 2500;
% Y = X.*cos(2*pi*(fo.*tX)+phi_b);
% 
% figure(1777);
% plot(tX,Y)
% 
% %B3
% 
% Sy = (var(Y)/T).*((abs(S)).^2);
% 
% figure(18);
% semilogy(Fx,Sy,'blue')
% set(gcf,'color', 'w');
% 
% hold on; 
% y1=Sy;
% y2=Sy;
% 
% SS=1/4*(y1+y2);
% 
% figure(19);
% semilogy(Fx,SS,'red')
% set(gcf,'color', 'w');
% 
% t2=-A*T+0:Ts:A*T+N*T-Ts;
% Ptotal=0;
% 
% for K=1:500
% theta = rand(1,1)*2*pi;
% X = conv(X_delta,phi)*Ts;
% X1= X.';
% Y=X.*cos(2*pi*fo*t2 + theta);
% S = fftshift(fft(Y, Nf)*Ts);
% Ttotal = length(Y)*Ts;
% Py = ((abs(S)).^2)/Ttotal;
% Ptotal = Ptotal+Py;
% end
% 
% Pfinal=Ptotal/K(end);
% 
% figure(20);
% plot(Fx,Pfinal,'blue');
% set(gcf,'color', 'w');
% 
% figure(21);
% semilogy(Fx,Pfinal,'blue');
% set(gcf,'color', 'w');















