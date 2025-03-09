clear; clc; close all;

%Α1 SRRC in time domain
T = 1e-3; over = 10; Ts = T/over; A = 4;
a_values = [0, 0.5, 1];

%cell array for the pulses
all_phi = cell(length(a_values),1);

figure('Name','Time-domain SRRC');
hold on; grid on;
for i=1:length(a_values)
    a = a_values(i);
    [phi, t] = srrc_pulse(T, over, A, a);
    plot(t, phi, 'DisplayName', sprintf('\\alpha=%.1f', a));
    
    % Save to use in A2
    all_phi{i} = phi;
end
xlabel('Time (s)'); ylabel('p(t)');
legend('Location','best');
title('SRRC pulses for different roll-offs');


%Α2
N_f = 1024; 
Fs = 1/Ts;
f = linspace(-Fs/2, Fs/2, N_f);
c = 1e-6;

figure('Name','Frequency-domain SRRC');

% (α) linear scale
subplot(2,1,1); hold on; grid on;
for i=1:length(a_values)
    phi = all_phi{i};
    P_f = fftshift( fft(phi, N_f) ) * Ts;
    plot(f, abs(P_f).^2, 'DisplayName', sprintf('\\alpha=%.1f', a_values(i)));
end
xlabel('f (Hz)'); ylabel('|Φ(f)|^2');
title('Energy spectral density (linear scale)');
legend('Location','best');

% (β) Log scale
subplot(2,1,2); hold on; grid on;
for i=1:length(a_values)
    phi = all_phi{i};
    P_f = fftshift( fft(phi, N_f) ) * Ts;
    semilogy(f, abs(P_f).^2, 'DisplayName', sprintf('\\alpha=%.1f', a_values(i)));
    yline(c, 'r--', 'Threshold c');
end
xlabel('f (Hz)'); ylabel('|Φ(f)|^2 (log scale)');
title('Energy spectral density (semilogy scale)');
legend('Location','best');
