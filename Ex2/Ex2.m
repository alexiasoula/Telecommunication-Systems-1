clear all; close all; clc

%A1
T= 10^(-2);
over = 10;
Ts = T/over;
A=4;
a=0.5;

[phi,t] = srrc_pulse(T,over,A,a);

Nf=2048;
Fs = 1/Ts;
f_axis = linspace(-Fs/2,Fs/2-Fs/Nf,Nf);
F_axis = Fs*f_axis;

PHI = fftshift(fft(phi,Nf))*Ts;
XF_abs = abs(PHI);
PHI_psd = XF_abs.^2;

figure(Name="A.1")
DrawSeminology(F_axis,PHI_psd, ...
    "Power Spectral Denstity",'Frequency','Amplitude');

%A.2
N = 100;
X = bits_to_2PAM(N);

figure(Name="A.2 bits to 2 PAM");
stem(X);
grid on;
title("Bits to 2PAM");
xlabel('Time');
ylabel('Amplitude');

X_delta = 1/Ts*upsample(X, over);
t_delta = [0 : Ts : (N * T) - Ts];

X_t = conv(X_delta,phi)*Ts;
t_conv = t(1)+t_delta(1):Ts:t(end)+t_delta(end);
figure(Name="A.2 Conv")
DrawPlot(t_conv,X_t,"Signal X(t)", ...
    "Time","Amplitude");

S_x = (var(X)/T)*PHI_psd;

%A.3
Px_F = calculateP(X_t,Nf,Ts,t_conv);

figure(Name="A.3 P_x(F)");
subplot(2,1,1);
DrawPlot(F_axis,Px_F,"Px(F) of signal X(t)",'Frequency','Amlitude');

subplot(2,1,2);
DrawSeminology(F_axis,Px_F,"Px(F) of signal X(t) [semilogy]", ...
    'Frequency','Amlitude');

k=500;

X_tests = zeros(k,Nf);
for i=1:k
    X_test = bits_to_2PAM(N);
    X_delta_test = 1/Ts*upsample(X_test, over);
    X_t_test = conv(X_delta_test,phi)*Ts;
    XF_psd_test = calculateP(X_t,Nf,Ts,1);
    X_tests(i,:)=XF_psd_test;
end

Sx_tests = mean(X_tests);

figure(Name="A.3 After tests");
semilogy(F_axis,S_x,'b');
hold on;
semilogy(F_axis,Sx_tests,'r');
hold off;
grid on;
title("Power Spectral Density");
xlabel('Frequency');
ylabel('Amplitude');
legend("Theoretical","Tests");

%A.4
X4 = bits_to_4PAM(N/2);
X4_delta = 1/Ts*upsample(X, over);
t4_delta = [0 : Ts : (N/2 * T) - Ts];
X4_t = conv(X_delta,phi)*Ts;
t4_conv = t(1)+t_delta(1):Ts:t(end)+t_delta(end);

figure(Name="A.4 Conv")
DrawPlot(t4_conv,X4_t,"Singal X(t)","Time","Amlitude");

Px_F4 = calculateP(X4_t,Nf,Ts,t4_conv);

figure(Name="A.4 Px(f)")
DrawSeminology(F_axis,Px_F4,"Px(F) of signal X(t)", ...
    'Frequency','Amplitude');



