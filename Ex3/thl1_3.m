clear all; close all; clc;
A = 1 ;
N     = 100;
T     = 1e-2;
over  = 10;
Ts    = T/over;
Fs    = 1/Ts;
beta  = 0.35;       % SRRC roll‑off factor
span  = 6;          % Filter length in symbols
roll = 0.5;
F0 = 200; %4
SNRdB = 15;   %7

bit_seq = generate_bits(N);    % 1. 4N random bits
X       = bits_to_PSK_16(bit_seq); % 2. 2×N 
XI      = X(1,:);          % {X_I,n}
XQ      = X(2,:);          % {X_Q,n}

% SRRC filter
[phi, t_phi] = srrc_pulse(T, over, A, roll);
% Upsampling
XI_up = (1/Ts)*upsample(XI, over);
XQ_up = (1/Ts)*upsample(XQ, over);

tx_delta_I = linspace(0, N*T, length(XI_up));
tx_delta_Q = linspace(0, N*T, length(XQ_up)); % same duration
% morfopiisi palmoy
sI = conv(XI_up, phi)* Ts; % I(t)
sQ = conv(XQ_up, phi)* Ts; % Q(t)

% time axes
t_conv = [tx_delta_I(1) + t_phi(1) : Ts : tx_delta_I(end) + t_phi(end) - Ts];

figure("Name","SRRC‑shaped Waveforms","NumberTitle","off");
subplot(2,1,1)
plot(t_conv, sI, 'b'); hold on;
plot(t_conv, sQ, 'r');
xlabel('Time (s)'); ylabel('Amplitude');
title('Baseband waveforms after SRRC shaping');
legend('s_I(t)','s_Q(t)'); grid on;

% periodiagramma
[pI, f] = periodogram(sI, [], 4096, Fs, 'centered');
[pQ, ~] = periodogram(sQ, [], 4096, Fs, 'centered');

subplot(2,1,2)
plot(f, 10*log10(pI), 'b'); hold on;
plot(f, 10*log10(pQ), 'r');
xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
title('Periodograms of I and Q branches');
legend('I‑branch','Q‑branch'); grid on;

%4
t = t_conv;
carrier_cos = cos(2*pi*F0*t);
carrier_sin = sin(2*pi*F0*t);

XI_t = sI .* carrier_cos;        % I(t) = s_I(t)·cos(2πF0 t)
XQ_t = -sQ .* carrier_sin;       % Q(t) = -s_Q(t)·sin(2πF0 t)

figure('Name','Passband Modulated Waveforms');
subplot(2,1,1)
plot(t, XI_t, 'b', t, XQ_t, 'r'); grid on;
legend('X_I(t)','X_Q(t)'); xlabel('Time (s)'); ylabel('Amplitude');
title(sprintf('Passband waveforms (F_0 = %d Hz)', F0));

[pXI,f2] = periodogram(XI_t, [], 4096, Fs, 'centered');
[pXQ,~]  = periodogram(XQ_t, [], 4096, Fs, 'centered');
subplot(2,1,2)
plot(f2, 10*log10(pXI), 'b', f2, 10*log10(pXQ), 'r'); grid on;
legend('X_I(t)','X_Q(t)'); xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
title('Periodograms');

%5
X_t = XI_t + XQ_t; %input signal

figure('Name','Channel Input');
subplot(2,1,1)
plot(t, X_t); grid on;
xlabel('Time (s)'); ylabel('Amplitude');
title('Channel input X(t)');

[pX,f3] = periodogram(X_t, [], 4096, Fs, 'centered');
subplot(2,1,2)
plot(f3, 10*log10(pX)); grid on;
xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
title('Periodogram of X(t)');

%7
% σ_W^2 = 1 / (Ts * 10^(SNRdB/10))
sigmaW2 = 1 / (Ts * 10^(SNRdB/10));

W_t = sqrt(sigmaW2) * randn(1, length(t));   % leukos Gaussian
Y_t = X_t + W_t;                             % output signal

figure('Name','Received signal Y(t)');
subplot(2,1,1)
plot(t, Y_t); grid on;
xlabel('Time (s)'); ylabel('Amplitude');
title(sprintf('Received signal Y(t)  (SNR = %d dB)', SNRdB));

[pY,f4] = periodogram(Y_t, [], 4096, Fs, 'centered');
subplot(2,1,2)
plot(f4, 10*log10(pY)); grid on;
xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
title('Periodogram of Y(t)');

%8
rxI = 2 * Y_t .* carrier_cos;
rxQ = -2 * Y_t .* carrier_sin;

figure('Name','Receiver Mixer Outputs');
subplot(2,1,1)
plot(t, rxI, 'b', t, rxQ, 'r'); grid on;
legend('r_I(t)','r_Q(t)'); xlabel('Time (s)'); ylabel('Amplitude');
title('Outputs after mixing down to baseband');

[prI,frI] = periodogram(rxI, [], 4096, Fs, 'centered');
[prQ,~]   = periodogram(rxQ, [], 4096, Fs, 'centered');
subplot(2,1,2)
plot(frI, 10*log10(prI), 'b', frI, 10*log10(prQ), 'r'); grid on;
legend('r_I','r_Q'); xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
title('Periodograms');

%9

yI = conv(rxI, phi)* Ts;
yQ = conv(rxQ, phi)* Ts;
t_conv_new = t(1) + t_phi(1) : Ts : t(end) + t_phi(end);
figure('Name','Receiver matched filter outputs');
subplot(2,1,1);
plot(t_conv_new, yI, 'b'); hold on; plot(t_conv_new, yQ, 'r');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;
title('Outputs after matched SRRC filter');
legend('y_I(t)','y_Q(t)');

subplot(2,1,2);
[PyI,fm] = periodogram(yI, [], 4096, Fs, 'centered');
[PyQ,~ ] = periodogram(yQ, [], 4096, Fs, 'centered');
plot(fm, 10*log10(PyI), 'b'); hold on;
plot(fm, 10*log10(PyQ), 'r'); grid on;
xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
title('Periodograms after matched filter');
legend('y_I','y_Q');

% A10 
A = 1;
idx_start = 2*A*over + 1;
idx_end = idx_start + (N-1) * over;

idx = idx_start : over : idx_end;
    

YI_samples = yI(idx);
YQ_samples = yQ(idx);

Y = [YI_samples; YQ_samples];

figure('Name','Scatterplot of received symbols');
scatter(YI_samples, YQ_samples, 'filled');
xlabel('In-phase'); ylabel('Quadrature');
title('Scatterplot after matched filter and sampling');
grid on;


%11
[X_est, bit_seq_est] = detect_PSK_16(Y);
num_errors_sumbols = symbol_errors(X_est,X);

%12
num_errors_bits = bit_errors(bit_seq_est,bit_seq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monte Carlo 1

N = 100;           
K = 1000;               
SNRdB = -2:2:24;        

numSNR = length(SNRdB);
num_symbol_errors = zeros(1,numSNR);
num_bit_errors = zeros(1,numSNR);

for idx = 1:numSNR
    snr = SNRdB(idx);
    for k = 1:K
        
        bit_seq = generate_bits(N);     
        X = bits_to_PSK_16(bit_seq);

        sigma2_W = 1/(Ts*10^(snr/10));
        sigma2_N = Ts*sigma2_W/2;
        noise = sqrt(sigma2_N)*randn(size(X));
        Y = X + noise;

        [est_X, est_bit_seq] = detect_PSK_16(Y);

        num_symbol_errors(idx) = num_symbol_errors(idx) + symbol_errors(est_X, X);
        num_bit_errors(idx) = num_bit_errors(idx) + bit_errors(est_bit_seq, bit_seq);
    end
end

Pe_symbol = num_symbol_errors/(K*N);
Pe_bit = num_bit_errors/(K*4*N);

figure;
semilogy(SNRdB, Pe_symbol, 'o-');
xlabel('SNR (dB)'); ylabel('Symbol Error Probability');
title('Monte Carlo εκτίμηση Pe(symbol) για 16-PSK');
grid on;


%2
M = 16;
Pe_symbol_bound = 2*qfunc(sqrt(2*10.^(SNRdB/10))*sin(pi/M));
hold on;
semilogy(SNRdB, Pe_symbol_bound, 'r--');
legend('Monte Carlo', 'Άνω φράγμα');

%3
Pe_bit_bound = Pe_symbol/log2(M);

figure;
semilogy(SNRdB, Pe_bit, 'o-'); hold on;
semilogy(SNRdB, Pe_bit_bound, 'r--');
xlabel('SNR (dB)'); ylabel('Bit Error Probability');
title('Monte Carlo εκτίμηση Pe(bit) για 16-PSK');
legend('Monte Carlo', 'Κάτω φράγμα');
grid on;