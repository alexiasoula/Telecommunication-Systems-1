clear all;
close all;
clc;

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
Nf = 4096;
A = 2;

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

t_delta = 0 : Ts : N*T-Ts;

% filtering with the SRRC pulse
XI_t = conv(XI_delta, phi) * Ts;
XQ_t = conv(XQ_delta, phi) * Ts;

% time of convolution
t_conv = t_delta(1) + t_phi(1) : Ts : t_delta(end) + t_phi(end);

% plots
figure
plot(t_conv, XI_t);
xlabel('sec');
ylabel('X_I(t)')
title('A3. X_I(t) after phi filter')

figure
plot(t_conv, XQ_t);
xlabel('sec');
ylabel('X_Q(t)')
title('A3. X_Q(t) after phi filter')

% periodograms
F = linspace(-Fs/2, Fs/2-Fs/Nf, Nf);
XI_f = fftshift( fft(XI_t, Nf) ) * Ts;
XQ_f = fftshift( fft(XQ_t, Nf) ) * Ts;
Tsig = t_conv(end) - t_conv(1);

P_I = abs(XI_f).^2 / Tsig;
P_Q = abs(XQ_f).^2 / Tsig;

% plot periodograms
figure;
plot(F, P_I);
xlabel('Freq (Hz)'); ylabel('P_{X_I}');
title('3) Periodogram of X_I');

figure;
plot(F, P_Q);
xlabel('Freq (Hz)'); ylabel('P_{X_Q}');
title('3) Periodogram of X_Q');


%% A4

Fo = 200;

XI_mod = 2 * XI_t .* cos(2*pi*Fo*t_conv);
XQ_mod = (-2) * XQ_t .* sin(2*pi*Fo*t_conv);


% plots
figure
plot(t_conv, XI_mod);
xlabel('sec');
ylabel('X_{Imod}(t)')
title('A4. X_I(t) after modulation')

figure
plot(t_conv, XQ_mod);
xlabel('sec');
ylabel('X_{Qmod}(t)')
title('A3. X_Q(t) after modulation')

% periodograms
F  = linspace(-Fs/2, Fs/2-Fs/Nf, Nf);
PI = abs(fftshift(fft(XI_mod,Nf))*Ts).^2 / (Tsig);
PQ = abs(fftshift(fft(XQ_mod,Nf))*Ts).^2 / (Tsig);


% plot periodograms
figure; plot(F, PI);
title('4) Periodogram of Modulated X_I'); xlabel('Freq (Hz)'); ylabel('P_I(f)');

figure; plot(F, PQ);
title('4) Periodogram of Modulated X_Q'); xlabel('Freq (Hz)'); ylabel('P_Q(f)');

%% A5

X_tx = XI_mod + XQ_mod; % transmitter

% plot
figure
plot(t_conv, X_tx);
xlabel('sec');
ylabel('X(t)');
title('A5. X transmitted')


X_tx_F = fftshift(abs(Ts*fft(X_tx,Nf)).^2);
P_x = X_tx_F / Tsig;

% plot
figure
plot(F,P_x);
xlabel('Freq(Hz)')
ylabel('P_X')
title('A5. X transmitted periodogram')
grid on

%% A6

Y = X_tx;

%% A7

% add gaussian noise at the output
SNR_dB   = 20;
sigma2_W = 1/(Ts * 10^(SNR_dB/10));
sigma2_N = Ts*sigma2_W/2;
W_t      = sqrt(sigma2_W)*randn(1, length(t_conv)); % white gaussian noise
Y     = X_tx + W_t; % noisy channel output

%% A8

YI_recv = Y .* cos(2*pi*Fo*t_conv);
YQ_recv = (-1) * Y .* sin(2*pi*Fo*t_conv);

% plots
figure
plot(t_conv, YI_recv)
xlabel('Time(s)')
ylabel('Y_{Irecv}')
title('A8. YI received and demodulated')
grid on

figure
plot(t_conv, YQ_recv)
xlabel('Time(s)')
ylabel('Y_{Qrecv}')
title('A8. YQ received and demodulated')
grid on

% periodograms

YI_F = fftshift(abs(Ts*fft(YI_recv,Nf)).^2);
YQ_F = fftshift(abs(Ts*fft(YQ_recv,Nf)).^2);

P_YI = YI_F/Tsig;
P_YQ = YQ_F/Tsig;

% plot periodograms
figure
plot(F,P_YI);
xlabel('Freq(Hz)')
ylabel('P_YI')
title('A8. periodogram YI')
grid on

figure
plot(F,P_YQ);
xlabel('Freq(Hz)')
ylabel('P_YQ')
title('A8. periodogram YQ')
grid on

%% A9

t_filter = t_conv(1) + t_phi(1) : Ts : t_conv(end) + t_phi(end);

Tsig = t_filter(end) - t_filter(1);

YI_phi = Ts*conv(YI_recv, phi);
YQ_phi = Ts*conv(YQ_recv, phi);

% plots
figure
plot(t_filter, YI_phi)
xlabel('sec')
title('A9. YI after phi filter')
grid on;

figure
plot(t_filter, YQ_phi)
xlabel('sec')
title('A9. YQ after phi filter')
grid on;

% periodograms
P_YI = fftshift(abs(Ts*fft(YI_phi,Nf)).^2)/Tsig;
P_YQ = fftshift(abs(Ts*fft(YQ_phi,Nf)).^2)/Tsig;

% plots
figure
plot(F,P_YI);
xlabel('Frequency(Hz)')
title('A9. periodogram YI_{phi}');
grid on

figure
plot(F,P_YQ);
xlabel('Frequency(Hz)')
title('A9. periodogram YQ_{phi}');
grid on

%% A10

YI_k = zeros(N,1);
YQ_k = zeros(N,1);

idx_y = 1;
for i = 2*A*over+1: over: length(t_filter)-2*A*over
    YI_k(idx_y) = YI_phi(i);
    YQ_k(idx_y) = YQ_phi(i);
    idx_y = idx_y + 1;
end

Yk = [YI_k, YQ_k];

scatterplot(Yk);
grid on;

%% A11, A12, A13

[est_X, est_bit_seq] = detect_PSK_16(Yk');

num_of_symbols_errors = symbol_errors(est_X, X)
num_of_bits_errors = bit_error(est_bit_seq, bit_seq)

%% B1
% Monte - Carlo

SNR_dB = -2:2:24;
K = 1000;
N = 100;
M = 16;
bps = log2(M);

num_of_symbols=N*K;
num_of_bits=4*N*K;


Pser_exp = zeros(1,length(SNR_dB));
Pber_exp = zeros(1,length(SNR_dB));

Pser_theory = zeros(1,length(SNR_dB));
Pber_theory = zeros(1,length(SNR_dB));

% parameters
T = 0.01;
over = 10;
Ts = T/over;
Fs = 1/Ts;
roll = 0.5; % roll-off factor
Nf = 2048;
A = 2;


for j = 1:length(SNR_dB)
    sum_sym_error = 0;
    sum_bit_error = 0;
    for i = 1:K
        bit_seq = (sign(randn(4*N,1)) + 1)/2; % 4N random bits (0 or 1)
        X = bits_to_PSK_16(bit_seq); % 2xN matrix
        XI = X(1,:);
        XQ = X(2,:);

        % SRRC pulse phi(t)
        [phi, t_phi] = srrc_pulse(T, over, A, roll);

        % delta-trains
        XI_delta = (1/Ts)*upsample(XI, over);
        XQ_delta = (1/Ts)*upsample(XQ, over);

        % define delta time axes

        t_delta = 0 : Ts : N*T-Ts;

        % filtering with the SRRC pulse
        XI_t = conv(XI_delta, phi) * Ts;
        XQ_t = conv(XQ_delta, phi) * Ts;

        % time of convolution
        t_conv = t_delta(1) + t_phi(1) : Ts : t_delta(end) + t_phi(end);

        Fo = 200;

        XI_mod = 2 * XI_t .* cos(2*pi*Fo*t_conv);
        XQ_mod = (-2) * XQ_t .* sin(2*pi*Fo*t_conv);

        X_tx = XI_mod + XQ_mod; % transmitter

        sigma2_W = 1/(Ts * 10^(SNR_dB(j)/10));
        sigma2_N = Ts*sigma2_W/2;
        W_t      = sqrt(sigma2_W)*randn(1, length(t_conv)); % white gaussian noise
        Y     = X_tx + W_t; % noisy channel output

        YI_recv = Y .* cos(2*pi*Fo*t_conv);
        YQ_recv = (-1) * Y .* sin(2*pi*Fo*t_conv);
        t_filter = t_conv(1) + t_phi(1) : Ts : t_conv(end) + t_phi(end);

        Tsig = t_filter(end) - t_filter(1);

        YI_phi = Ts*conv(YI_recv, phi);
        YQ_phi = Ts*conv(YQ_recv, phi);

        YI_k = zeros(N,1);
        YQ_k = zeros(N,1);

        idx_y = 1;
        for i = 2*A*over+1: over: length(t_filter)-2*A*over
            YI_k(idx_y) = YI_phi(i);
            YQ_k(idx_y) = YQ_phi(i);
            idx_y = idx_y + 1;
        end

        Yk = [YI_k, YQ_k];
        [est_X, est_bit_seq] = detect_PSK_16(Yk');

        num_of_symbols_errors = symbol_errors(est_X, X);
        num_of_bits_errors = bit_error(est_bit_seq, bit_seq);

        sum_sym_error = sum_sym_error + num_of_symbols_errors;
        sum_bit_error = sum_bit_error + num_of_bits_errors;
    end
    Pser_exp(j) = sum_sym_error / num_of_symbols;
    Pber_exp(j) = sum_bit_error / num_of_bits;
    
    SNR = 10^(SNR_dB(j)/10);

    Pser_theory(j) = 2*Q(sqrt(2*(SNR))*sin(pi/16));
    Pber_theory(j) = Pser_theory(j)/4;
end

%% B2

figure
semilogy(SNR_dB, Pser_exp);
hold on;
semilogy(SNR_dB, Pser_theory, 'r');
xlabel('SNR (db)')
ylabel('log(probabilty)')
title('B2. probability of symbol error')

figure
semilogy(SNR_dB, Pber_exp);
hold on;
semilogy(SNR_dB, Pber_theory, 'r');
xlabel('SNR (db)')
ylabel('log(probabilty)')
title('B2. probability of bit error')

 



