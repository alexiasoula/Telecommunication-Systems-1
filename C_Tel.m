clear; clc; close all;
% function for C.2 1
function [X] = bits_to_2PAM(b)
    for i = 1 : length(b)
        if b(i) == 0
            X(i) = 1;
        else 
            X(i) = -1;
        end
    end
end
%C.1
N = 100;

T = 10^(-3);
over = 10;
a = 0.5;
A = 4;

b = (sign(randn(N, 1)) + 1)/2;


%C.2 1 Call the function
X = bits_to_2PAM(b);

% C.2 2 
T_s = T / over;

X_delta = 1/T_s * upsample(X, over);

t_d = [0: T_s : (N * T) - T_s];

figure('Name','C.2 2 Draw X_delta');
stem(t_d, X_delta);
title('Delta pulse');
xlabel('Time');
ylabel('Amplitude');
grid on;

% C.2 3
[phi, t] = srrc_pulse(T, over, A, a);

X_t = conv(X_delta,phi)*T_s;
t_conv = [min(t_d) + min(t) : T_s : max(t_d) + max(t)];

figure('Name','C.2 3 Draw X_t');
plot(t_conv, X_t);
title('Conv Delta pulse with phi');
xlabel('Time');
ylabel('Amplitude');
grid on;

% C.2 4
% Create φ(-t)
phi_negative = phi(numel(phi):-1:1);
t_negative =  t(numel(t):-1:1);

Z = conv(X_t, phi_negative)*T_s;
t_z = [min(t_negative) + min(t_conv) : T_s : max(t_negative) + max(t_conv)];

kT = [0: N - 1] * T;

figure('Name','C.2 4 Draw Z_t');
plot(t_z, Z);
hold on;
title('Z(t)');
xlabel('Time');
ylabel('Amplitude');
grid on;
stem(kT, X);
hold off;







