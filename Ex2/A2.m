function [X, X_delta,t_delta, X_t, t_conv, S_x] = A2(N,over,T,Ts,phi,t,PHI_psd)
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

end

