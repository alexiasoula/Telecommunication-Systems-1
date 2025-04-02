function [phi,t, PHI_psd, F_axis] = A1(T, over, A, a, Fs, Nf, Ts)
    [phi,t] = srrc_pulse(T,over,A,a);
    
    
    f_axis = linspace(-Fs/2,Fs/2-Fs/Nf,Nf);
    F_axis = Fs*f_axis;
    
    PHI_psd = calculateP(phi,Nf,Ts,1);

    figure(Name="A.1")
    DrawSeminology(F_axis,PHI_psd, ...
        "Power Spectral Denstity",'Frequency','Amplitude');

end

