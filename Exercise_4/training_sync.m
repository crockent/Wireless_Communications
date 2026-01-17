function [C_d, corr_Y, d_opt_corr] = training_sync(Z_r, training_s, over, Ntr, analog_channel)
    B_sym = 3;  % SRRC support ~2*B symbols [web:1 from simulation]
    d_total = length(analog_channel) - B_sym * over;
    need_len_tr = d_total + (Ntr-1)*over;
    if length(Z_r) < need_len_tr
        Z_r = [Z_r; zeros(need_len_tr - length(Z_r), 1)];
    end
    C_d = zeros(d_total, 1);
    corr_Y = zeros(d_total, 1);
    for dd = 1:d_total
        idx = dd:over: (dd + (Ntr-1)*over);
        C_d(dd) = sum(conj(training_s) .* Z_r(idx));
        corr_Y(dd) = abs(C_d(dd));
    end
    [~, d_opt_corr] = max(corr_Y);
end
