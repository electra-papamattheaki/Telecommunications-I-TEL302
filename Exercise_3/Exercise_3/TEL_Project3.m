clear all;
close all;
clc;

% 1.
N = 100;
bit_seq = (sign(randn(4*N, 1)) + 1)/2;

% 2. 
X = bits_to_PSK_16(bit_seq);

Xi = X(1,:);
Xq = X(2,:);

% Plot PSK Symbols
figure()
plot(Xi,Xq,'o')
%plot(Xi,Xq)
set(gcf,'color', 'w');
grid on;
title('PSK Plot')

% 3.
% Parameters
T = 0.01;
over = 10;
Ts = T/over;
a = 1;
A = 1; 
Nf = 4096;
Fs = 1/Ts;
Fx = -Fs/2:Fs/Nf:Fs/2-Fs/Nf;

% upsampling Xi, Xq
X_delta_i = 1/Ts * upsample(Xi, over);
X_delta_q = 1/Ts * upsample(Xq, over);

t = 0:Ts:N*Ts*over-Ts;
l = t(1):Ts:t(end);

figure();
stem(t,X_delta_i);
title('Waveform of X_i,n');
grid on;
set(gcf,'color', 'w');

figure();
stem(t,X_delta_q);
title('Waveform of X_q,n');
grid on;
set(gcf,'color', 'w');

% creation of SRRC Pulse
[phi, t1] = srrc_pulse(T, over, A, a);

% Convolution
X_I = conv(X_delta_i,phi)*Ts; 
tX_I = t1(1)+t(1):Ts:t1(end)+t(end);

X_Q = conv(X_delta_q,phi)*Ts;
tX_Q = t1(1)+t(1):Ts:t1(end)+t(end);

figure();
plot(tX_I,X_I);
title('Output Waveform of X_i(t)');
ylabel('X_i(t)');
xlabel('Time(s)');
grid on; 
set(gcf,'color', 'w');

figure();
plot(tX_Q,X_Q);
title('Output Waveform of X_q(t)');
ylabel('X_q(t)');
xlabel('Time(s)');
set(gcf,'color', 'w');

% Plot the periodograms
Ttotal_I = tX_I(end) - tX_I(1);
S_I = fftshift(fft(X_I,Nf)* Ts);
P_I = ((abs(S_I)) .^ 2) / Ttotal_I;

Ttotal_Q = tX_Q(end) - tX_Q(1);
S_Q = fftshift(fft(X_Q,Nf)* Ts);
P_Q = ((abs(S_Q)) .^ 2) / Ttotal_Q;

figure();
plot(Fx, P_I , 'blue');
title('Periodogram of X_I using plot');
ylabel('Spectrum');
xlabel('Frequency(Hz)');
grid on; 
set(gcf,'color', 'w');

figure();
semilogy(Fx, P_I , 'blue');
title('Periodogram of X_I using semilogy');
ylabel('Spectrum');
xlabel('Frequency(Hz)');
grid on; 
set(gcf,'color', 'w');

figure();
plot(Fx, P_Q , 'blue');
title('Periodogram of X_Q using plot');
ylabel('Spectrum');
xlabel('Frequency(Hz)');
grid on; 
set(gcf,'color', 'w');

figure();
semilogy(Fx, P_Q , 'blue');
title('Periodogram of X_Q using semilogy');
ylabel('Spectrum');
xlabel('Frequency(Hz)');
grid on; 
set(gcf,'color', 'w');

% 4. 

fo = 200;
X_I_mod = 2 * X_I .* cos(2*pi*(fo.*tX_I));
X_Q_mod = (-2) * X_Q .* sin(2*pi*(fo.*tX_Q));

X_I_2 = X_I.*X_I_mod;
X_Q_2 = X_Q.*X_Q_mod;

figure();
%plot(tX_I,X_I_2);
plot(tX_I,X_I_mod);
title('Modulated X_i(t) Waveform')
xlabel('Time(s)');
ylabel('Amplitude');
grid on;
set(gcf,'color', 'w');

figure()
%plot(tX_Q,X_Q_2);
plot(tX_Q,X_Q_mod);
title('Modulated X_q(t) Waveform')
xlabel('time(s)');
ylabel('Amplitude');
grid on;
set(gcf,'color', 'w');

% Plot the periodograms
S_I_mod = fftshift(fft(X_I_mod,Nf)* Ts);
P_I_mod = ((abs(S_I_mod)) .^ 2) / Ttotal_I;

S_Q_mod = fftshift(fft(X_Q_mod,Nf)* Ts);
P_Q_mod = ((abs(S_Q_mod)) .^ 2) / Ttotal_Q;

figure();
plot(Fx, P_I_mod , 'blue');
title('Periodogram of modulate X_i using plot');
ylabel('Spectrum');
xlabel('Frequency(Hz)');
grid on; 
set(gcf,'color', 'w');

figure();
semilogy(Fx, P_I_mod , 'blue');
title('Periodogram of modulate X_i using semilogy');
ylabel('Spectrum');
xlabel('Frequency(Hz)');
grid on;
set(gcf,'color', 'w');

figure();
plot(Fx, P_Q_mod , 'blue');
title('Periodogram of modulate X_q using plot');
ylabel('Spectrum');
xlabel('Frequency(Hz)');
grid on;
set(gcf,'color', 'w');

figure();
semilogy(Fx, P_Q_mod , 'blue');
title('Periodogram of modulate X_q using semilogy');
ylabel('Spectrum');
xlabel('Frequency(Hz)');
grid on;
set(gcf,'color', 'w');

% 5. 
X_mod = X_I_mod + X_Q_mod;
figure();
plot(tX_I,X_mod);
title('Waveform of Modulated Signal X');
xlabel('time(s)');
ylabel('Amplitude');
grid on; 
set(gcf,'color', 'w');

S_X_mod = fftshift(fft(X_mod,Nf)* Ts);
P_X_mod = ((abs(S_X_mod)) .^ 2) / Ttotal_I;

figure();
plot(Fx, P_X_mod , 'blue');
title('Periodogram of Modulated Signal X using plot');
xlabel('Time(s)');
ylabel('Amplitude');
grid on;
set(gcf,'color', 'w');

figure();
semilogy(Fx, P_X_mod , 'blue');
title('Periodogram of Modulated Signal X using semilogy');
xlabel('Time(s)');
ylabel('Amplitude');
grid on;
set(gcf,'color', 'w');

% 6. 
% Considered ideal!

% 7.
% Gaussian noise
SNR_db = 20;
var_W = (10*(A^2))/ (Ts* 10^(SNR_db/10)); 
Wt = sqrt(var_W).*randn(1,length(X_mod)); 
% Adding noise
Y = X_mod + Wt;

figure();
subplot(2,1,1)
plot(tX_I,X_mod);
title('Signal X without Noise');
xlabel('time(s)');
ylabel('Amplitude');
subplot(2,1,2)
plot(tX_I,Y);
title('Signal Y with Noise');
xlabel('time(s)');
ylabel('Amplitude');
grid on; 
set(gcf,'color', 'w');

% 8.
X_Wi_mod = Y .* cos(2*pi*(fo.*tX_I));
X_Wq_mod = (-1) * Y .* sin(2*pi*(fo.*tX_Q));

figure();
plot(tX_Q,X_Wi_mod);
title('Demodulated Signal X_i');
xlabel('time(s)');
ylabel('Amplitude');
grid on; 
set(gcf,'color', 'w');

figure()
plot(tX_Q,X_Wq_mod);
title('Demodulated Signal X_q');
xlabel('time(s)');
ylabel('Amplitude');
set(gcf,'color', 'w');

% Plot the periodograms
Ttotal_Wi = tX_I(end) - tX_I(1);
S_Wi = fftshift(fft(X_Wi_mod,Nf)* Ts);
P_Wi = ((abs(S_Wi)) .^ 2) / Ttotal_Wi;

Ttotal_Wq = tX_Q(end) - tX_Q(1);
S_Wq = fftshift(fft(X_Wq_mod,Nf)* Ts);
P_Wq = ((abs(S_Wq)) .^ 2) / Ttotal_Wq;

figure();
plot(Fx, P_Wi , 'blue');
title('Periodogram of Demodulated X_i using plot');
ylabel('Spectrum');
xlabel('Frequency(Hz)');
grid on;
set(gcf,'color', 'w');

figure();
semilogy(Fx, P_Wi , 'blue');
title('Periodogram of Demodulated X_i using semilogy');
ylabel('Spectrum');
xlabel('Frequency(Hz)');
grid on;
set(gcf,'color', 'w');

figure();
plot(Fx, P_Wq , 'blue');
title('Periodogram of Demodulated X_q using plot');
ylabel('Spectrum');
xlabel('Frequency(Hz)');
grid on;
set(gcf,'color', 'w');

figure();
semilogy(Fx, P_Wq , 'blue');
title('Periodogram of Demodulated X_q using semilogy');
ylabel('Spectrum');
xlabel('Frequency(Hz)');
grid on;
set(gcf,'color', 'w');

% 9. 
X_Wi_conv = conv(X_Wi_mod,phi)*Ts; 
tX_I1 = tX_I(1)+t1(1):Ts:tX_I(end)+t1(end);

X_Wq_conv = conv(X_Wq_mod,phi)*Ts;
tX_Q1 = tX_Q(1)+t1(1):Ts:tX_Q(end)+t1(end);

figure();
plot(tX_I1,X_Wi_conv);
title('Filtered Signal X_i');
xlabel('Time(s)');
ylabel('Amplitude');
grid on; 
set(gcf,'color', 'w');

figure();
plot(tX_Q1,X_Wq_conv);
title('Filtered Signal X_q');
xlabel('Time(s)');
ylabel('Amplitude');
grid on; 
set(gcf,'color', 'w');

Ttotal_Wi_conv = tX_I1(end) - tX_I1(1);
S_Wi_conv = fftshift(fft(X_Wi_conv,Nf)* Ts);
P_Wi_conv = ((abs(S_Wi_conv)) .^ 2) / Ttotal_Wi_conv;

Ttotal_Wq_conv = tX_Q1(end) - tX_Q1(1);
S_Wq_conv = fftshift(fft(X_Wq_conv,Nf)* Ts);
P_Wq_conv = ((abs(S_Wq_conv)) .^ 2) / Ttotal_Wq_conv;

figure();
plot(Fx, P_Wi_conv , 'blue');
title('Periodogram of Filtered X_i using plot');
ylabel('Spectrum');
xlabel('Frequency(Hz)');
grid on; 
set(gcf,'color', 'w');

figure();
semilogy(Fx, P_Wi_conv , 'blue');
title('Periodogram of Filtered X_i using semilogy');
ylabel('Spectrum');
xlabel('Frequency(Hz)');
grid on;
set(gcf,'color', 'w');

figure();
plot(Fx, P_Wq_conv , 'blue');
title('Periodogram of Filtered X_q using plot');
ylabel('Spectrum');
xlabel('Frequency(Hz)');
grid on;
set(gcf,'color', 'w');

figure();
semilogy(Fx, P_Wq_conv , 'blue');
title('Periodogram of Filtered X_i using semilogy');
ylabel('Spectrum');
xlabel('Frequency(Hz)');
grid on;
set(gcf,'color', 'w');

% 10. 
YI = X_Wi_conv(2*over+1 :over: (2+N)*over);
YQ = X_Wq_conv(2*over+1 :over: (2+N)*over);

Y_sampled = YI+1i*YQ;

scatterplot(Y_sampled)
set(gcf,'color', 'w');

% 11. 
[Y_est, est_bit_seq] = detect_PSK_16(Y_sampled);

% 12.
num_of_symbol_errors = symbol_errors(Y_est, X);

% 13. 
num_of_bit_errors = bit_errors(est_bit_seq, bit_seq); 

%% Part 2. 

% 1. 

K = 200;
SNR_db_B = -2:2:24;
symbolErrB = 0;
bitsErrB = 0;
symbolErrB_Total = 0;
bitsErrB_Total = 0;
total_symbols = N*K;
total_bits = 4*N*K;
P_error_symbol =  ones(1,length(SNR_db_B));
P_error_bit = ones(1,length(SNR_db_B));

for i = 1:length(SNR_db_B)
    for j = 1:K
        [symbolErrB, bitsErrB] = PSK_16(SNR_db_B(i));
        symbolErrB_Total = symbolErrB_Total + symbolErrB;
        bitsErrB_Total = bitsErrB_Total + bitsErrB;        
    end

    P_error_symbol(i) = symbolErrB_Total / total_symbols;
    P_error_bit(i) = bitsErrB_Total / total_bits;
    symbolErrB_Total = 0;
    bitsErrB_Total = 0;
end

% Theoritical calculation
sw = 10.*(A^2)./(Ts .* 10.^(SNR_db_B ./ 10));
sn =(Ts*sw)/2;
Theor_symbolErr = 3*Q(A./ sqrt(sn));
Theor_bitErr = Theor_symbolErr/log2(16);

% 2. 
figure()
semilogy(SNR_db_B, P_error_symbol,'blue');
hold on;
semilogy(SNR_db_B, Theor_symbolErr,'red');
legend('Experimental','Theoretical');
xlabel('SNR in dB');
ylabel('Symbol Error Rate for 16-PSK');
legend('Experimental','Theoretical')
grid on; 
set(gcf,'color', 'w');

% 3.
figure()
semilogy(SNR_db_B, P_error_bit,'blue');
hold on;
semilogy(SNR_db_B, Theor_bitErr,'red');
xlabel('SNR in dB');
ylabel('Bit Error Rate for 16-PSK');
legend('Experimental','Theoretical')
grid on;
set(gcf,'color', 'w');



