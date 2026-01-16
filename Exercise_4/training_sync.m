function [C_d, corr_Y, d_opt_corr] = training_sync(Z_r, trainig_s_known, over, Ntr,analog_channel)
% TRAINING_SYNC_METRIC  Computes training-based synchronization metric |C_d|
%   [C_d, corr_Y, d_opt_corr] = training_sync_metric(Z_r, trainig_s_known, over, Ntr)
%   - Z_r: matched filtered received signal (complex)
%   - trainig_s_known: known training symbols (complex, length Ntr)
%   - over: oversampling factor
%   - Ntr: number of training symbols
%   Returns:
%   - C_d: complex correlation values for each delay (length d_total)
%   - corr_Y: |C_d| metric
%   - d_opt_corr: optimal delay argmax |C_d|

A_est = 3;
d_total = length(analog_channel) - A_est * over;  % max possible delay

% Ensure Z_r is long enough for max delay
need_len_tr = d_total + (Ntr-1)*over;
if length(Z_r) < need_len_tr
    Z_r = [Z_r; zeros(need_len_tr - length(Z_r), 1)];
end

C_d = zeros(d_total, 1);
corr_Y = zeros(d_total, 1);
for dd = 1:d_total
    idx = dd : over : dd + (Ntr-1)*over;
    C_d(dd) = sum(conj(trainig_s_known) .* Z_r(idx));
    corr_Y(dd) = abs(C_d(dd));
end

[~, d_opt_corr] = max(corr_Y);  % optimal synchronization delay

end
