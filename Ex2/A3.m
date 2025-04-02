function A3(a,X_t,Nf,T,t_conv,F_axis,N,over,phi,S_x)
    Ts = T / over; 

    Px_F = calculateP(X_t,Nf,Ts,t_conv);
    
    figure(Name="A."+a+" P_x(F)");
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
        XF_psd_test = calculateP(X_t_test,Nf,Ts,1);
        X_tests(i,:)=XF_psd_test;
    end
    
    Sx_tests = mean(X_tests);
    
    figure(Name="A."+a+" After tests");
    semilogy(F_axis,S_x,'b');
    hold on;
    semilogy(F_axis,Sx_tests,'r');
    hold off;
    grid on;
    title("Power Spectral Density");
    xlabel('Frequency');
    ylabel('Amplitude');
    legend("Theoretical","Tests");
end

