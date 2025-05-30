clear all;
close all;
clc;
%rng(2)

%% A1
N = 100;
bit_seq = (sign(randn(4*N,1)) + 1)/2; % 4N random bits (0 or 1)

%% A3

% parameters
T = 0.01;
over = 10;
Ts = T/over;
Fs = 1/Ts;
roll = 0.5; % roll-off factor
Nf = 2048;
A = 1;

% map the bit_seq to 16-PSK (Gray coded)
X = bits_to_PSK_16(bit_seq); % 2xN matrix
XI = X(1,:);
XQ = X(2,:);

% SRRC pulse phi(t)
[phi, t_phi] = srrc_pulse(T, over, A, roll);

% delta-trains
XI_delta = (1/Ts)*upsample(XI, over);
XQ_delta = (1/Ts)*upsample(XQ, over);

% define delta time axes

tx_delta_I = linspace(0, N*T, length(XI_delta));
tx_delta_Q = linspace(0, N*T, length(XQ_delta)); % same duration

% filtering with the SRRC pulse
XI_t = conv(XI_delta, phi) * Ts;
XQ_t = conv(XQ_delta, phi) * Ts;

%allagi edw
%step = (t_phi(end) + tx_delta_I(end) - (t_phi(1) + tx_delta_I(1))) / length(XI_t);
t_conv = [tx_delta_I(1) + t_phi(1) : Ts : tx_delta_I(end) + t_phi(end) - Ts];

% time-domain plots
figure;
plot(t_conv, XI_t);
xlabel('Time (s)'); ylabel('X_I(t)');
title('3) X_I(t) signal after SRRC');

figure;
plot(t_conv, XQ_t);
xlabel('Time (s)'); ylabel('X_Q(t)');
title('3) X_Q(t) signal after SRRC');

% periodograms
F = linspace(-Fs/2, Fs/2-Fs/Nf, Nf);
XI_f = fftshift( fft(XI_t, Nf) ) * Ts;
XQ_f = fftshift( fft(XQ_t, Nf) ) * Ts;
Tsig = t_conv(end) - t_conv(1);

P_I = abs(XI_f).^2 / Tsig;
P_Q = abs(XQ_f).^2 / Tsig;

figure;
plot(F, P_I);
xlabel('Freq (Hz)'); ylabel('P_{X_I}');
title('3) Periodogram of X_I');

figure;
semilogy(F, P_I);
xlabel('Freq (Hz)'); ylabel('P_{X_I} (log)');
title('3) Periodogram of X_I (log)');

figure;
plot(F, P_Q);
xlabel('Freq (Hz)'); ylabel('P_{X_Q}');
title('3) Periodogram of X_Q');

figure;
semilogy(F, P_Q);
xlabel('Freq (Hz)'); ylabel('P_{X_Q} (log)');
title('3) Periodogram of X_Q (log)');

%% A4

F0 = 200; % carrier frequency
t = t_conv; % time axis

XI_afterCos = 2 * XI_t .* cos(2*pi*F0*t);
XQ_afterSin = -2 * XQ_t .* sin(2*pi*F0*t);

% plot modulated XI,XQ
figure; plot(t, XI_afterCos);
title('4) Modulated X_I'); xlabel('Time (s)'); ylabel('X_I(t)·cos(2ðF_0t)');

figure; plot(t, XQ_afterSin);
title('4) Modulated X_Q'); xlabel('Time (s)'); ylabel('-X_Q(t)·sin(2ðF_0t)');

% periodograms
F  = linspace(-Fs/2, Fs/2-Fs/Nf, Nf);
PI = abs(fftshift(fft(XI_afterCos,Nf))*Ts).^2 / (t(end)-t(1));
PQ = abs(fftshift(fft(XQ_afterSin,Nf))*Ts).^2 / (t(end)-t(1));

figure; plot(F, PI);
title('4) Periodogram of Modulated X_I'); xlabel('Freq (Hz)'); ylabel('P_I(f)');

figure; plot(F, PQ);
title('4) Periodogram of Modulated X_Q'); xlabel('Freq (Hz)'); ylabel('P_Q(f)');

%% A5

X_t = XI_afterCos + XQ_afterSin;

% plot X(t)
figure;
plot(t_conv, X_t);
xlabel('Time (s)');
ylabel('X(t)');
title('5) Channel Input Signal X(t)');

% compute periodogram of X(t)
X_f  = fftshift( fft(X_t, Nf) ) * Ts;
Tsig = t_conv(end) - t_conv(1);
P_X  = abs(X_f).^2 / Tsig;

% plot the periodogram
figure;
plot(F, P_X);
xlabel('Frequency (Hz)');
ylabel('P_{X}');
title('5) Periodogram of X(t)');

% log-scale version
figure;
semilogy(F, P_X);
xlabel('Frequency (Hz)');
ylabel('P_{X} [log]');
title('5) Periodogram of X(t) (log scale)');

%% A6

% ideal channel so the output equals the input
Y_t = X_t; % channel output before noise

%% A7

% add white Gaussian noise at the output
SNR_dB   = 20;
sigma2_W = 1/(Ts * 10^(SNR_dB/10));
W_t      = sqrt(sigma2_W)*randn(size(X_t)); % white Gaussian noise
Y_t      = X_t + W_t; % noisy channel output


%% A8

YI_afterCos =  Y_t .* cos(2*pi*F0*t);
YQ_afterSin = -(Y_t .* sin(2*pi*F0*t));

% plot them
figure;
plot(t, YI_afterCos);
xlabel('Time (s)');
ylabel('Amplitude');
title('8) I-branch Output (before filtering)');

figure;
plot(t, YQ_afterSin);
xlabel('Time (s)');
ylabel('Amplitude');
title('8) Q-branch Output (before filtering)');

% plot their periodograms
P_YI = abs(fftshift(fft(YI_afterCos, Nf)) * Ts).^2 / Tsig;
P_YQ = abs(fftshift(fft(YQ_afterSin, Nf)) * Ts).^2 / Tsig;

figure;
plot(F, P_YI);
xlabel('Frequency (Hz)');
ylabel('Power');
title('8) Periodgram of I-branch Output');

figure;
plot(F, P_YQ);
xlabel('Frequency (Hz)');
ylabel('Power');
title('8) Periodgram of Q-branch Output');

%% A9

% apply the SRRC filter
YI_filtered = conv(YI_afterCos, phi) * Ts;
YQ_filtered = conv(YQ_afterSin, phi) * Ts;

% t_start = t_conv(1) + t_phi(1);  % Use t_conv from A3
% t_end = t_conv(end) + t_phi(end);
% t_conv_new  = linspace(t_start, t_end, length(YI_filtered));
t_conv_new = t(1) + t_phi(1) : Ts : t(end) + t_phi(end);

% plot them
figure;
plot(t_conv_new, YI_filtered);
xlabel('Time (s)');
ylabel('Y_{I,filtered}(t)');
title('9) I-branch After SRRC Filter');

figure;
plot(t_conv_new, YQ_filtered);
xlabel('Time (s)');
ylabel('Y_{Q,filtered}(t)');
title('9) Q-branch After SRRC Filter');

% plot the periodograms
Tsig_new = t_conv_new(end) - t_conv_new(1);

P_YI_f = abs(fftshift(fft(YI_filtered, Nf)) * Ts).^2 / Tsig_new;
P_YQ_f = abs(fftshift(fft(YQ_filtered, Nf)) * Ts).^2 / Tsig_new;

figure;
plot(F, P_YI_f);
xlabel('Frequency (Hz)');
ylabel('Power');
title('9) Periodogram of I-branch After Filtering');

figure;
plot(F, P_YQ_f);
xlabel('Frequency (Hz)');
ylabel('Power');
title('9) Periodogram of Q-branch After Filtering');

%% A10  –
A = 1;
idx_start = 2*A*over + 1;
idx_end = idx_start + (N-1) * over;

idx = idx_start : over : idx_end;
YI_symbols = YI_filtered(idx);
YQ_symbols = YQ_filtered(idx);

Y = [YI_symbols; YQ_symbols];

figure;
scatter(YI_symbols, YQ_symbols, 'filled');
xlabel('In-Phase (I)');
ylabel('Quadrature (Q)');
title('Scatter Plot of Received Symbols');
grid on;
axis equal;

% ideal symbols (no noise)
M = 16;
theta = 2*pi*(0:M-1)/M;
I_ref = cos(theta);
Q_ref = sin(theta);

hold on;
scatter(I_ref, Q_ref, 80, 'r', 'o', 'LineWidth', 1.5); % red
legend('Received symbols', 'Ideal constellation');


%% A11 , A12, A13

[X_est, bit_seq_est] = detect_PSK_16(Y);

num_of_symbols_error = symbol_errors(X_est, X)
num_of_bits_error = bit_error(bit_seq_est, bit_seq)


%% B1
% Monte - Carlo

SNR_dB = -2:2:24;
K = 1000;
N = 200;
M = 16;
bps = log2(M);

num_of_symbols=N*K;
num_of_bits=4*N*K;

Pser_exp = zeros(1,length(SNR_dB));
Pber_expr = zeros(1,length(SNR_dB));

Pser_theory = zeros(1,length(SNR_dB));
Pber_theory = zeros(1,length(SNR_dB));

bit_seq = (sign(randn(4*N,1)) + 1)/2; % 4N random bits (0 or 1)

% parameters
T = 0.01;
over = 10;
Ts = T/over;
Fs = 1/Ts;
roll = 0.5; % roll-off factor
Nf = 2048;
A = 1;

sum_sym_error = 0;
sum_bit_error = 0;
for j =SNR_dB
for i=1:K
    % map the bit_seq to 16-PSK (Gray coded)
    X = bits_to_PSK_16(bit_seq); % 2xN matrix
    XI = X(1,:);
    XQ = X(2,:);

    % SRRC pulse phi(t)
    [phi, t_phi] = srrc_pulse(T, over, A, roll);

    % delta-trains
    XI_delta = (1/Ts)*upsample(XI, over);
    XQ_delta = (1/Ts)*upsample(XQ, over);

    % define delta time axes

    tx_delta_I = linspace(0, N*T, length(XI_delta));
    tx_delta_Q = linspace(0, N*T, length(XQ_delta)); % same duration

    % filtering with the SRRC pulse
    XI_t = conv(XI_delta, phi) * Ts;
    XQ_t = conv(XQ_delta, phi) * Ts;

    t_conv = tx_delta_I(1) + t_phi(1) : Ts : tx_delta_I(end) + t_phi(end) - Ts;

    F0 = 200; % carrier frequency
    t = t_conv; % time axis

    XI_afterCos = 2 * XI_t .* cos(2*pi*F0*t);
    XQ_afterSin = -2 * XQ_t .* sin(2*pi*F0*t);


    X_t = XI_afterCos + XQ_afterSin;

    sigma2_W = 1/(Ts * 10^(i/10));
    W_t      = sqrt(sigma2_W)*randn(size(X_t)); % white Gaussian noise
    Y_t      = X_t + W_t; % noisy channel output

    YI_afterCos =  Y_t .* cos(2*pi*F0*t);
    YQ_afterSin = -(Y_t .* sin(2*pi*F0*t));

    YI_filtered = conv(YI_afterCos, phi) * Ts;
    YQ_filtered = conv(YQ_afterSin, phi) * Ts;

    % t_start = t_conv(1) + t_phi(1);  % Use t_conv from A3
    % t_end = t_conv(end) + t_phi(end);
    % t_conv_new  = linspace(t_start, t_end, length(YI_filtered));
    t_conv_new = t(1) + t_phi(1) : Ts : t(end) + t_phi(end);

    idx_start = 2*A*over + 1;
    idx_end = idx_start + (N-1) * over;

    idx = idx_start : over : idx_end;
    YI_symbols = YI_filtered(idx);
    YQ_symbols = YQ_filtered(idx);

    Y = [YI_symbols; YQ_symbols];

    [X_est, bit_seq_est] = detect_PSK_16(Y);

    num_of_symbols_error = symbol_errors(X_est, X);
    num_of_bits_error = bit_error(bit_seq_est, bit_seq);

    sum_sym_error = sum_sym_error + num_of_symbols_error;
    sum_bit_error = sum_bit_error + num_of_bits_error;

end
end

Pser_exp = sum_sym_error / num_of_symbols
Pber_expr = sum_bit_error / num_of_bits

%% B2

