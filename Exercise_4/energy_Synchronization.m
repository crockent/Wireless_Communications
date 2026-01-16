function [Ed,pos_max_energy] = energy_Synchronization(Y_matched, over, N_symbols,analog_channel)
    % ENERGY_SYNC_STATISTIC: Στατιστικό ενέργειας για symbol timing sync (logic only)
    % Inputs:
    %   Y_matched: matched filter output 
    %   over: oversampling factor (Ts/T=10)
    %   N_symbols: αριθμός συμβόλων για μέσο ενέργειας 
    % Outputs:
    %   pos_max_energy: d_opt = argmax E_d (1-based index)
    %   energy_stat: vector E_d για d=1:d_total

    A_est = 3;  % approx SRRC support B=3 -> 6T, adjust based on phi
    d_total = length(analog_channel) - A_est * over;


    need_len = d_total + (N_symbols-1)*over;
    num_trailing = max(0, need_len - length(Y_matched));
    Y_pad = [Y_matched; complex(zeros(num_trailing,1))];  % complex zero pad

    Ed = zeros(d_total, 1);
    for dd = 1:d_total
        idx = dd : over : dd + (N_symbols-1)*over;
        Ed(dd) = sum(abs(Y_pad(idx)).^2);
    end

    [~, pos_max_energy] = max(Ed);
    
end

