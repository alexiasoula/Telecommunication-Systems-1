clear; clc; close all;

%parameters
T=1e-3;
A=10;
over = 10;
a_values=[0,0.5,1];
k_values=[0,1,2,3];

for a = a_values
    [phi_base, t_base]= srrc_pulse(T, over, A, a);
    Ts = T/over;

    figure('Name',sprintf('SRRC alpha=%.2f',a));

    for ki= 1:length(k_values)
        k=k_values(ki);
        
        %How many samples are in time kT
        sample_shift= round(k*T/Ts);
        
        %if sample_shift outside bounds, care for the shift
        phi_shifted = zeros(size(phi_base));
        if sample_shift >= 0
            % right shift: first sample_shift samples are zeros
            phi_shifted(sample_shift+1:end) = phi_base(1:end-sample_shift);
        else
            % left shift
            phi_shifted(1:end+sample_shift) = phi_base(-sample_shift+1:end);
        end
         % 1. φ(t), φ(t-kT) on the same plot
        subplot(2, length(k_values), ki)
        plot(t_base, phi_base, 'b','LineWidth',1.2); hold on; grid on;
        plot(t_base, phi_shifted, 'r--','LineWidth',1.2);
        xlabel('t (sec)');
        ylabel('\phi(t)');
        legend('\phi(t)','\phi(t - kT)','Location','best');
        title(sprintf('\\alpha=%.1f, k=%d', a, k));
        
        % 2. φ(t)*φ(t-kT)
        product = phi_base .* phi_shifted;
        subplot(2, length(k_values), length(k_values)+ki)
        plot(t_base, product, 'k','LineWidth',1.2); grid on;
        xlabel('t (sec)'); ylabel('\phi(t)\phi(t-kT)');
        title('Product');
        
    end
    % 3. calculate integrals k=0,1,2,3
    fprintf('\n=== alpha = %.2f ===\n', a);
    for k = k_values
        sample_shift = round(k*T / Ts);
        phi_shifted = zeros(size(phi_base));
        if sample_shift >= 0
            phi_shifted(sample_shift+1:end) = phi_base(1:end-sample_shift);
        else
            phi_shifted(1:end+sample_shift) = phi_base(-sample_shift+1:end);
        end
        
        % Product
        product = phi_base .* phi_shifted;
        
        % Integral (trapz)
        integral_val = trapz(t_base, product);
        
        fprintf('   k=%d -> Integral ~ %g\n', k, integral_val);
    end
end