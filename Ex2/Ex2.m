clear all; close all; clc

%A1
T= 10^(-2);
over = 10;
Ts = T/over; 
A=4;
a=0.5;

Nf=2048;
Fs = 1/Ts;

[phi, t, PHI_psd, F_axis] = A1(T,over,A,a,Fs,Nf);

%A.2
N = 100;
[X,X_delta,t_delta, X_t, t_conv, S_x] = A2(N,over,T,phi,t,PHI_psd);

%A.3

A3("3",X_t,Nf,T,t_conv,F_axis,N,over,phi,S_x);

%A.4
X4 = bits_to_4PAM(N/2);
X4_delta = 1/Ts*upsample(X, over);
t4_delta = [0 : Ts : N/2 * T - Ts];
X4_t = conv(X_delta,phi)*Ts;
t4_conv = linspace(t(1)+t4_delta(1),t(end)+t4_delta(end),length(X4_t));
%t4_conv = t(1)+t4_delta(1):Ts:t(end)+t4_delta(end);

figure(Name="A.4 Conv")
DrawPlot(t4_conv,X4_t,"Singal X(t)","Time","Amlitude");

Px_F4 = calculateP(X4_t,Nf,Ts,t4_conv);

figure(Name="A.4 Px(f)")
DrawSeminology(F_axis,Px_F4,"Px(F) of signal X(t)", ...
    'Frequency','Amplitude');

k=500;

X_tests4 = zeros(k,Nf);
for i=1:k
    X_test4 = bits_to_4PAM(N/2);
    X_delta_test4 = 1/Ts*upsample(X_test4, over);
    X_t_test4 = conv(X_delta_test4,phi)*Ts;
    XF_psd_test4 = calculateP(X_t_test4,Nf,Ts,1);
    X_tests4(i,:)=XF_psd_test4;
end

Sx_tests4 = mean(X_tests4);
S_x4 = (var(X4)/T)*PHI_psd;

figure(Name="A.4 Theoretical - Tests");
semilogy(F_axis,S_x4);
hold on;
DrawSeminology(F_axis,Sx_tests4,"Power Spectral Density", ...
    "Frequency","Amlitude")
hold off;
legend("Theoretical","Tests");

%A5
A3("5",X_t,Nf,2*T,t_conv,F_axis,N,2*over,phi,S_x);
