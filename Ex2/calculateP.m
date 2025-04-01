function [Px_f] = calculateP(X,Nf,Ts,t_conv)
    X_F = fftshift(fft(X,Nf))*Ts;

    T_total = max(t_conv)- min(t_conv)+1;

    Px_f = (abs(X_F)).^2 / T_total;
end