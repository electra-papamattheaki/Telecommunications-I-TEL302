function [num_of_symbol_errors , num_of_bit_errors] = PSK_16(SNR_db)

% 1.
N = 100;
bit_seq = (sign(randn(4*N, 1)) + 1)/2;

% 2. 
X = bits_to_PSK_16(bit_seq);

% 3.
Xi = X(1,:);
Xq = X(2,:);

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

% creation of SRRC Pulse
[phi, t1] = srrc_pulse(T, over, A, a);

% Convolution
X_I = conv(X_delta_i,phi)*Ts; 
tX_I = t1(1)+t(1):Ts:t1(end)+t(end);

X_Q = conv(X_delta_q,phi)*Ts;
tX_Q = t1(1)+t(1):Ts:t1(end)+t(end);

% Plot the periodograms
Ttotal_I = tX_I(end) - tX_I(1);
S_I = fftshift(fft(X_I,Nf)* Ts);
P_I = ((abs(S_I)) .^ 2) / Ttotal_I;

Ttotal_Q = tX_Q(end) - tX_Q(1);
S_Q = fftshift(fft(X_Q,Nf)* Ts);
P_Q = ((abs(S_Q)) .^ 2) / Ttotal_Q;

% 4. 

fo = 200;
X_I_mod = 2 * X_I .* cos(2*pi*(fo.*tX_I));
X_Q_mod = (-2) * X_Q .* sin(2*pi*(fo.*tX_Q));

X_I_2 = X_I.*X_I_mod;
X_Q_2 = X_Q.*X_Q_mod;



% Plot the periodograms
S_I_mod = fftshift(fft(X_I_mod,Nf)* Ts);
P_I_mod = ((abs(S_I_mod)) .^ 2) / Ttotal_I;

S_Q_mod = fftshift(fft(X_Q_mod,Nf)* Ts);
P_Q_mod = ((abs(S_Q_mod)) .^ 2) / Ttotal_Q;

% 5. 
X_mod = X_I_mod + X_Q_mod;

S_X_mod = fftshift(fft(X_mod,Nf)* Ts);
P_X_mod = ((abs(S_X_mod)) .^ 2) / Ttotal_I;

% 6. 
% Considered ideal!

% 7.
% Gaussian noise
var_W = (10*(A^2))/ (Ts* 10^(SNR_db/10)); 
Wt = sqrt(var_W).*randn(1,length(X_mod)); 
% Adding noise
Y = X_mod + Wt;

% 8.
X_Wi_mod = Y .* cos(2*pi*(fo.*tX_I));
X_Wq_mod = (-1) * Y .* sin(2*pi*(fo.*tX_Q));

Ttotal_Wi = tX_I(end) - tX_I(1);
S_Wi = fftshift(fft(X_Wi_mod,Nf)* Ts);
P_Wi = ((abs(S_Wi)) .^ 2) / Ttotal_Wi;

Ttotal_Wq = tX_Q(end) - tX_Q(1);
S_Wq = fftshift(fft(X_Wq_mod,Nf)* Ts);
P_Wq = ((abs(S_Wq)) .^ 2) / Ttotal_Wq;

% 9. 
X_Wi_conv = conv(X_Wi_mod,phi)*Ts; 
tX_I1 = tX_I(1)+t1(1):Ts:tX_I(end)+t1(end);

X_Wq_conv = conv(X_Wq_mod,phi)*Ts;
tX_Q1 = tX_Q(1)+t1(1):Ts:tX_Q(end)+t1(end);

Ttotal_Wi_conv = tX_I1(end) - tX_I1(1);
S_Wi_conv = fftshift(fft(X_Wi_conv,Nf)* Ts);
P_Wi_conv = ((abs(S_Wi_conv)) .^ 2) / Ttotal_Wi_conv;

Ttotal_Wq_conv = tX_Q1(end) - tX_Q1(1);
S_Wq_conv = fftshift(fft(X_Wq_conv,Nf)* Ts);
P_Wq_conv = ((abs(S_Wq_conv)) .^ 2) / Ttotal_Wq_conv;

% 10. 
YI = X_Wi_conv(2*over+1 :over: (2+N)*over);
YQ = X_Wq_conv(2*over+1 :over: (2+N)*over);

Y_sampled = YI+1i*YQ;

% 11. 

[Y_est, est_bit_seq] = detect_PSK_16(Y_sampled);

% 12.
num_of_symbol_errors = symbol_errors(Y_est, Y_sampled);

% 13. 
num_of_bit_errors = bit_errors(est_bit_seq, bit_seq); 

end