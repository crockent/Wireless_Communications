clear all; clc; close all;

%% =========================
%% PART 1: 
%% =========================

%% 1.1)
T = 1;          % symbol period
over = 10;      % oversampling factor 
Ts = T/over;    % sample period
roll_off = 0.4; % roll-off
A = 3;          % half-length of SRRC in symbol periods

% Construct SRRC pulse
g(:,1) = srrc_pulse2(T, Ts, A, roll_off);

M = 10;      % length of discrete equivalent channel
sA2 = 2;     % 4-QAM power

% Generate data
N = 200;       % # of symbols per packet
N_tr = 50;     % # of training symbpls per packet

% Training symbols
train_bits = (sign(randn(2*N_tr,1)) + 1)/2;
data_bits  = (sign(randn(2*(N-N_tr),1)) + 1)/2;

train_data = bits_to_4QAM(train_bits); % size N_tr x 1
data       = bits_to_4QAM(data_bits);  % size (N-N_tr) x 1

% Generation of packet
packet = [train_data; data];

% Shaping filter

% Pass the upampled packet from the SRRC filter 
% - you may also plot the real and the imaginary parts of X_Tx
packet_up = upsample(packet, over);
X_Tx = conv(packet_up, g) ;

% figure; 
% plot([real(X_Tx) imag(X_Tx)]), title('Real/Imag part of X_Tx') 

%% 1.2) - 1.8) 
for case_id = 1:2

    if case_id == 1
        c = 1;  % ideal channel: delta(t)
        case_name = 'Ideal channel';
    else
        c = (randn + 1i*randn)/sqrt(2);  % flat fading: c*delta(t), c in CN(0,1)
        case_name = sprintf('Flat fading: c=%.2f%+.2fj', real(c), imag(c));
    end

    % 1.2 / 1.8) Output of physical channel
    Y_Rx = c * X_Tx;
    %figure; 
    %plot([real(Y_Rx) imag(Y_Rx)]), title(['Real/Imag part of Y\_Rx - ' case_name])

    % 1.3) Matched filter
    Y_Rx_matched = conv(Y_Rx, g) * Ts;
    %figure; 
    %plot([real(Y_Rx_matched) imag(Y_Rx_matched)]), title(['Real/Imag part of Y\_Rx\_matched - ' case_name])

    % 1.4) Composite channel
    h_composite = c * conv(g, g) * Ts;

    figure;
    plot([real(h_composite) imag(h_composite)]), title(['Real/Imag part of composite channel - ' case_name])

    figure; plot(abs(h_composite));
    grid on; title(['|h_{composite}[n]| - ' case_name]);
    xlabel('sample index n'); ylabel('|h|');

    % 1.5) Synchronization via output energy
    d_total = length(h_composite) - A*over;

    need_len = d_total + (N-1)*over;
    num_of_trailing_zeros = max(0, need_len - length(Y_Rx_matched));
    Y_Rx_matched = [Y_Rx_matched; zeros(num_of_trailing_zeros,1)];

    energy_Y = zeros(d_total,1);
    for dd = 1:d_total
        idx = dd : over : (dd + (N-1)*over);
        energy_Y(dd,1) = sum(abs(Y_Rx_matched(idx)).^2);
    end
    [~,pos_max_energy_Y] = max(energy_Y);

    figure; hold on;
    plot([real(h_composite) imag(h_composite)])
    sample_idx = pos_max_energy_Y:over:pos_max_energy_Y+(M-1)*over;
    sample_idx = sample_idx(sample_idx <= length(h_composite));   % safety cap
    stem(sample_idx, real(h_composite(sample_idx)), '--b')
    stem(sample_idx, imag(h_composite(sample_idx)), '--r')
    title(['Composite channel & sampled taps (Part 1.5) - ' case_name]); grid on;

    % 1.6) Synchronization via training (C_d)
    need_len_tr = d_total + (N_tr-1)*over;
    if length(Y_Rx_matched) < need_len_tr
        Y_Rx_matched = [Y_Rx_matched; zeros(need_len_tr - length(Y_Rx_matched),1)];
    end

    C_d    = zeros(d_total,1);
    corr_Y = zeros(d_total,1);
    for dd = 1:d_total
        idx = dd : over : (dd + (N_tr-1)*over);
        C_d(dd) = sum(conj(train_data) .* Y_Rx_matched(idx));
        corr_Y(dd) = abs(C_d(dd));
    end

    figure;
    plot(corr_Y); grid on;
    title(['Training-based metric |C_d| - ' case_name]); xlabel('d (MATLAB index)'); ylabel('|C_d|');

    h_hat = C_d / (N_tr * sA2);

    E_d = zeros(d_total,1);
    for dd = 1:d_total
        sample_idx = dd : over : (dd + (M-1)*over);
        sample_idx = sample_idx(sample_idx <= length(h_hat));  % keep valid
        E_d(dd) = sum(abs(h_hat(sample_idx)).^2);
    end

    figure;
    plot(E_d); grid on;
    title(['E_d from h_hat - ' case_name]); xlabel('d (MATLAB index)'); ylabel('E_d');

    Lplot = min(d_total, length(h_composite));
    figure; hold on; grid on;
    plot(corr_Y(1:Lplot)/max(corr_Y(1:Lplot)+eps), 'LineWidth', 1.2);
    plot(abs(h_composite(1:Lplot))/max(abs(h_composite(1:Lplot))+eps), '--', 'LineWidth', 1.2);
    legend('|C_d| (norm)','|h_{composite}| (norm)');
    title(['Training metric vs |h_{composite}| (norm) - ' case_name]);
    xlabel('d (MATLAB index)'); ylabel('normalized');

    % 1.7) Symbol-spaced output sequence (length N)
    pos_sync_E = pos_max_energy_Y;

    [~,pos_max_Cd] = max(corr_Y);
    pos_sync_T = pos_max_Cd;

    need_len_Y = max(pos_sync_E, pos_sync_T) + (N-1)*over;
    if length(Y_Rx_matched) < need_len_Y
        Y_Rx_matched = [Y_Rx_matched; zeros(need_len_Y - length(Y_Rx_matched),1)];
    end

    Y_E = Y_Rx_matched(pos_sync_E : over : pos_sync_E + (N-1)*over);
    Y_T = Y_Rx_matched(pos_sync_T : over : pos_sync_T + (N-1)*over);

    figure;
    plot(real(Y_E), imag(Y_E), 'o'); axis equal; grid on;
    title(['Unequalized data (timing from energy, 1.5) - ' case_name]);
    xlabel('Re'); ylabel('Im');

    figure;
    plot(real(Y_T), imag(Y_T), 'o'); axis equal; grid on;
    title(['Unequalized data (timing from training, 1.6) - ' case_name]);
    xlabel('Re'); ylabel('Im');

    disp(['[ ' case_name ' ]  pos_energy, pos_Cd, diff = ' ...
          num2str(pos_max_energy_Y) '  ' num2str(pos_max_Cd) '  ' num2str(pos_max_energy_Y-pos_max_Cd)]);

end

%% =========================
%% PART 2: 
%% =========================

Dfs = 1e-3;                  % CFO at sampling rate: ÄFs = ÄF*Ts = 10^(-3)
phi_list = [0, 2*pi*rand];   % case1: phi=0, case2: random phi (same random value within this run)
phi_tag  = {'phi=0', 'phi random'};

for kphi = 1:2

    phi0 = phi_list(kphi);

    % Ideal channel: y[n] = x[n]
    c = 1;
    Y_Rx = c * X_Tx;

    % Apply CFO at sample-rate: y_CFO[n] = exp(j(2ðÄFs n + ö)) * y[n]
    n_samp   = (0:length(Y_Rx)-1).';
    CFO_term = exp(1j*(2*pi*Dfs*n_samp + phi0));
    Y_Rx_CFO = Y_Rx .* CFO_term;

    % Matched filter
    Y_Rx_matched = conv(Y_Rx_CFO, g) * Ts;

    % Composite channel 
    h_composite = c * conv(g, g) * Ts;

    % Synchronization via output energy
    d_total = length(h_composite) - A*over;

    need_len = d_total + (N-1)*over;
    num_of_trailing_zeros = max(0, need_len - length(Y_Rx_matched));
    Y_Rx_matched = [Y_Rx_matched; zeros(num_of_trailing_zeros,1)];

    energy_Y = zeros(d_total,1);
    for dd = 1:d_total
        idx = dd : over : (dd + (N-1)*over);
        energy_Y(dd) = sum(abs(Y_Rx_matched(idx)).^2);
    end
    [~, pos_sync] = max(energy_Y);

    % Symbol-spaced output (length N)
    Y = Y_Rx_matched(pos_sync : over : pos_sync + (N-1)*over);

    figure;
    plot(real(Y), imag(Y), 'o'); axis equal; grid on
    title(sprintf('Part 2: Unequalized symbol-spaced output (with CFO, %s)', phi_tag{kphi}));
    xlabel('Re'); ylabel('Im');

    % CFO estimation using training symbols 
    train_matr = toeplitz(train_data(M:N_tr), train_data(M:-1:1).');
    B = train_matr * pinv(train_matr' * train_matr) * train_matr';

    v_axis = -0.03 : 1e-5 : 0.03;
    f = zeros(length(v_axis),1);

    y_seg = Y(M:N_tr);
    Lseg  = length(y_seg);
    l_vec = (0:Lseg-1).';

    for i_v = 1:length(v_axis)
        exp_i_v = exp(1j * 2*pi * v_axis(i_v) * l_vec);
        G_i_v   = diag(exp_i_v);
        f(i_v)  = real( y_seg' * G_i_v * B * G_i_v' * y_seg );
    end

    figure;
    plot(v_axis, f), grid on
    title(sprintf('Part 2: Statistic for CFO estimation (training-based, %s)', phi_tag{kphi}));
    xlabel('v (cycles/symbol)'); ylabel('statistic');

    [~, pos_v] = max(f);
    v_est = v_axis(pos_v);

    v_true = over * Dfs;
    fprintf('[Part 2 - %s] v_true = over*Dfs = %.6f,  v_est = %.6f,  phi = %.3f rad\n', ...
            phi_tag{kphi}, v_true, v_est, phi0);

    % CFO correction at symbol level 
    n_sym = (0:length(Y)-1).';
    Z = Y .* exp(-1j * 2*pi * v_est * n_sym);

    figure;
    plot(real(Z), imag(Z), 'o'); axis equal; grid on
    title(sprintf('Part 2: Symbol-spaced output AFTER CFO correction (%s)', phi_tag{kphi}));
    xlabel('Re'); ylabel('Im');

    % 1-tap LS channel estimate 
    Z_tr = Z(1:N_tr);
    h_ls = (train_data(1:N_tr)' * Z_tr) / (train_data(1:N_tr)' * train_data(1:N_tr));
    fprintf('[Part 2 - %s] LS channel estimate (1 tap): h_ls = %.4f%+.4fj, |h_ls|=%.4f\n', ...
            phi_tag{kphi}, real(h_ls), imag(h_ls), abs(h_ls));

    % Equalize 
    Z_eq = Z / h_ls;

    figure;
    plot(real(Z_eq), imag(Z_eq), 'o'); axis equal; grid on
    title(sprintf('Part 2: After CFO correction + channel equalization (%s)', phi_tag{kphi}));
    xlabel('Re'); ylabel('Im');

    % prints
    h0 = h_composite(2*A*over + 1);
    fprintf('[Part 2 - %s] h0 = %.4f%+.4fj\n', phi_tag{kphi}, real(h0), imag(h0));

end


%% =========================
%% PART 3:
%% =========================

% 3.1) Physical channel: 2-tap 
K = randi([1, 4*over]);
c0 = (randn + 1i*randn)/sqrt(2);
rho = 0.6;
c1 = rho*(randn + 1i*randn)/sqrt(2);

h_physical = zeros(K+1,1);
h_physical(1)   = c0;
h_physical(K+1) = c1;

fprintf('[Part 3.1] K=%d samples, |c0|=%.3f, |c1|=%.3f (rho=%.2f)\n', ...
        K, abs(c0), abs(c1), rho);

% 3.2–3.3) Output + matched filter + composite channel 
Y_Rx = conv(X_Tx, h_physical);
Y_Rx_matched = conv(Y_Rx, g) * Ts;

h_composite = conv(h_physical, conv(g, g)) * Ts;

figure; plot([real(h_composite) imag(h_composite)]), grid on
title(sprintf('Part 3: Real/Imag of composite channel (K=%d samples)', K));

figure; plot(abs(h_composite)), grid on
title('Part 3: |h_{composite}[n]|'); xlabel('n'); ylabel('|h|');

[hc_max, nhc] = max(abs(h_composite));
fprintf('[Part 3.3] len(Y_Rx_matched)=%d, len(h_composite)=%d, max|h|=%.3f at n=%d\n', ...
        length(Y_Rx_matched), length(h_composite), hc_max, nhc-1);

% 3.4) Training-correlation statistic corrd(d)=|C_d|
d_total = 4*A*over + 1 + K; 

need_len_tr = d_total + (N_tr-1)*over;
if length(Y_Rx_matched) < need_len_tr
    Y_Rx_matched = [Y_Rx_matched; zeros(need_len_tr - length(Y_Rx_matched),1)];
end

C_d   = zeros(d_total,1);
corrd = zeros(d_total,1);

for dd = 1:d_total
    idx = dd : over : (dd + (N_tr-1)*over);
    C_d(dd)   = sum(conj(train_data) .* Y_Rx_matched(idx));
    corrd(dd) = abs(C_d(dd));
end

[~, dC_idx] = max(corrd);
dC = dC_idx - 1; % 0-based

fprintf('[Part 3.4] peak |C_d| at dC=%d (0-based)\n', dC);

figure; plot(corrd); grid on
title('Part 3: corrd(d)=|C_d|'); xlabel('dd (MATLAB index)'); ylabel('|C_d|');

% 3.5) Compute E_d, find sync, and take symbol-spaced sequence (length N+M-1)

h_hat = C_d / (N_tr * sA2);

E_d = zeros(d_total,1);

% pad h_hat so dd+(M-1)*over is always valid
h_hat_pad = [h_hat; zeros((M-1)*over,1)];

for dd = 1:d_total
    sample_idx = dd : over : (dd + (M-1)*over);
    E_d(dd) = sum(abs(h_hat_pad(sample_idx)).^2);
end

figure; plot(E_d); grid on
title('Part 3.5: E_d from h_hat'); xlabel('dd (MATLAB index)'); ylabel('E_d');

% Synchronization
[~, pos_sync] = max(E_d);              
d_star = pos_sync - 1; % 0-based delay (samples)

fprintf('[Part 3.5] pos_sync=%d (MATLAB), d_star=%d (0-based), d_total=%d\n', ...
        pos_sync, d_star, d_total);
fprintf('[Check] dC_idx=%d, length(E_d)=%d\n', dC_idx, length(E_d));

% Symbol-spaced output sequence of length N+M-1 
need_len_Y = pos_sync + (N+M-2)*over;
if length(Y_Rx_matched) < need_len_Y
    Y_Rx_matched = [Y_Rx_matched; zeros(need_len_Y - length(Y_Rx_matched),1)];
end

Y = Y_Rx_matched(pos_sync : over : pos_sync + (N+M-2)*over); % length N+M-1

fprintf('[Part 3.5] len(Y)=%d (expected %d)\n', length(Y), N+M-1);

figure;
plot(real(Y), imag(Y), 'o'); axis equal; grid on
title('Part 3.5: Unequalized symbol-spaced output (length N+M-1)');
xlabel('Re'); ylabel('Im');

fprintf('[Part 3.5] phase = mod(pos_sync-1,over) = %d\n', mod(pos_sync-1, over));

% 3.5b) LS estimate of discrete equiv. channel (length M)

if N_tr < M
    error('3.5b: Need N_tr >= M to estimate an M-tap channel (LS).');
end

% Training convolution matrix 
train_matr = toeplitz(train_data(M:N_tr), train_data(M:-1:1).');   

% Observed output
y_ls = Y(M:N_tr);

% LS estimate
h_d_ls = (train_matr' * train_matr) \ (train_matr' * y_ls);

fprintf('[Part 3.5b] LS M-tap channel estimate computed: size=%dx1\n', length(h_d_ls));

% test
resid = y_ls - train_matr * h_d_ls;
fprintf('[Part 3.5b] LS residual energy = %.4e (mean |resid|^2 = %.4e)\n', ...
        sum(abs(resid).^2), mean(abs(resid).^2));

% Compare with "true" discrete equivalent from composite channel sampling 
need_len_hc = pos_sync + (M-1)*over;
if length(h_composite) < need_len_hc
    h_comp_pad = [h_composite; zeros(need_len_hc - length(h_composite), 1)];
else
    h_comp_pad = h_composite;
end

h_d_true = h_comp_pad(pos_sync : over : pos_sync + (M-1)*over);
h_d_true_plot = h_d_true;   

figure;
plot(0:M-1, abs(h_d_true_plot), 'o-'); grid on; hold on;
plot(0:M-1, abs(h_d_ls),   's-'); 
title('Part 3.5b: |h_d| true (from composite) vs LS estimate');
xlabel('tap index m'); ylabel('magnitude');
legend('|h_d true|', '|h_d LS|');

% PART 3.6)
h_d = h_d_ls(:); % same as LS channel estimate 
Mch = length(h_d);

L = 5*Mch;               
delta = Mch;              

H = toeplitz([h_d; zeros(L-1,1)], [h_d(1); zeros(L-1,1)]);

% Desired overall response: e_delta, with 1 at index (delta+1)
e_delta = zeros(Mch+L-1,1);
e_delta(delta+1) = 1;

% ZF (LS) solution
f_zf = (H' * H) \ (H' * e_delta);

fprintf('[Part 3.6] ZF equalizer computed: L=%d, delta=%d (0-based)\n', L, delta);

% check: overall impulse response
g_overall = conv(h_d, f_zf);   % length M+L-1

[pk, ipk] = max(abs(g_overall));
fprintf('[Part 3.6] max|conv(h_d,f)|=%.3f at n=%d (0-based), target delta=%d\n', ...
        pk, ipk-1, delta);

figure; stem(0:length(g_overall)-1, abs(g_overall)); grid on
title('Part 3.6: |conv(h_d, f_{ZF})| (should peak at n=\delta)');
xlabel('n'); ylabel('|.|');


% 3.7) 

% Filter symbol-spaced sequence with ZF equalizer
Z = conv(Y, f_zf); % length (N+M-1) + L - 1

% Discard first delta samples at equalizer output
start_idx = delta + 1; % MATLAB index 
end_idx   = start_idx + N - 1;

% Padding
if length(Z) < end_idx
    Z = [Z; zeros(end_idx - length(Z),1)];
end

Z_keep = Z(start_idx:end_idx);

fprintf('[Part 3.7] len(Z)=%d, taking Z_keep(%d:%d) => len=%d (want N=%d)\n', ...
        length(Z), start_idx, end_idx, length(Z_keep), N);

figure;
plot(real(Z_keep), imag(Z_keep), 'o');
axis equal; grid on
title('Part 3.7: After ZF equalization (discard \delta, keep N)');
xlabel('Re'); ylabel('Im');

% check
fprintf('[Part 3.7] mean |Z_keep|^2 = %.4f\n', mean(abs(Z_keep).^2));

%% =========================
%% PART 4: 
%% =========================

% % 4.1) Physical channel: 2-tap + CFO
% K = randi([1, 4*over]);
% c0 = (randn + 1i*randn)/sqrt(2);
% rho = 0.6;
% c1 = rho*(randn + 1i*randn)/sqrt(2);
% 
% h_physical = zeros(K+1,1);
% h_physical(1)   = c0;
% h_physical(K+1) = c1;

% CFO at sampling-rate: 
Dfs4 = 1e-4*(rand - 0.5);
phi4 = 2*pi*rand;

fprintf('[Part 4.1] K=%d, |c0|=%.3f, |c1|=%.3f, Dfs=%.3g, phi=%.3f rad\n', ...
        K, abs(c0), abs(c1), Dfs4, phi4);

% 4.2) Output 
Y_Rx = conv(X_Tx, h_physical);

% 4.3) Apply CFO at sample-rate
n_samp4 = (0:length(Y_Rx)-1).';
CFO_term4 = exp(1j*(2*pi*Dfs4*n_samp4 + phi4));
Y_Rx_CFO = Y_Rx .* CFO_term4;

% 4.4) Matched filter
Y_Rx_matched = conv(Y_Rx_CFO, g) * Ts;

% 4.5) Composite channel 
h_composite = conv(h_physical, conv(g,g)) * Ts;

figure; plot(abs(h_composite)), grid on
title('Part 4: |h_{composite}[n]|'); xlabel('n'); ylabel('|h|');

% 4.6) Training correlation corrd(d)=|C_d|  
d_total = 4*A*over + 1 + K;

need_len_tr = d_total + (N_tr-1)*over;
if length(Y_Rx_matched) < need_len_tr
    Y_Rx_matched = [Y_Rx_matched; zeros(need_len_tr - length(Y_Rx_matched),1)];
end

C_d   = zeros(d_total,1);
corrd = zeros(d_total,1);

for dd = 1:d_total
    idx = dd : over : (dd + (N_tr-1)*over);
    C_d(dd)   = sum(conj(train_data) .* Y_Rx_matched(idx));
    corrd(dd) = abs(C_d(dd));
end

[~, dC_idx] = max(corrd);
fprintf('[Part 4.6] peak |C_d| at dC=%d (0-based)\n', dC_idx-1);

figure; plot(corrd); grid on
title('Part 4: corrd(d)=|C_d| (with CFO)'); xlabel('dd (MATLAB index)'); ylabel('|C_d|');

% 4.7) Compute E_d from h_hat, find pos_sync, take symbol-spaced Y (length N+M-1)
h_hat = C_d / (N_tr * sA2);

E_d = zeros(d_total,1);
h_hat_pad = [h_hat; zeros((M-1)*over,1)]; 

for dd = 1:d_total
    sample_idx = dd : over : (dd + (M-1)*over);
    E_d(dd) = sum(abs(h_hat_pad(sample_idx)).^2);
end

[~, pos_sync] = max(E_d);
fprintf('[Part 4.7] pos_sync=%d (MATLAB), d_star=%d (0-based), phase=%d\n', ...
        pos_sync, pos_sync-1, mod(pos_sync-1,over));

figure; plot(E_d); grid on
title('Part 4: E_d from h_hat (with CFO)'); xlabel('dd (MATLAB index)'); ylabel('E_d');

need_len_Y = pos_sync + (N+M-2)*over;
if length(Y_Rx_matched) < need_len_Y
    Y_Rx_matched = [Y_Rx_matched; zeros(need_len_Y - length(Y_Rx_matched),1)];
end

Y = Y_Rx_matched(pos_sync : over : pos_sync + (N+M-2)*over);  % length N+M-1

figure; plot(real(Y), imag(Y), 'o'); axis equal; grid on
title('Part 4: Unequalized symbol-spaced output (with CFO)');
xlabel('Re'); ylabel('Im');

% 4.8) CFO estimation at symbol-rate 
% v_true = over*Dfs4 (cycles/symbol)
train_matr = toeplitz(train_data(M:N_tr), train_data(M:-1:1).');
B = train_matr * pinv(train_matr' * train_matr) * train_matr';

v_axis = -0.03 : 1e-5 : 0.03;
fstat = zeros(length(v_axis),1);

y_seg = Y(M:N_tr);
Lseg  = length(y_seg);
l_vec = (0:Lseg-1).';

for i_v = 1:length(v_axis)
    exp_i_v = exp(1j * 2*pi * v_axis(i_v) * l_vec);
    G_i_v   = diag(exp_i_v);
    fstat(i_v)  = real( y_seg' * G_i_v * B * G_i_v' * y_seg );
end

[~, pos_v] = max(fstat);
v_est = v_axis(pos_v);

fprintf('[Part 4.8] v_true=%.6f, v_est=%.6f\n', over*Dfs4, v_est);

figure; plot(v_axis, fstat), grid on
title('Part 4: CFO statistic (training-based)'); xlabel('v (cycles/symbol)'); ylabel('stat');

% 4.9) CFO correction at symbol-level
n_sym = (0:length(Y)-1).';
Zcfo = Y .* exp(-1j * 2*pi * v_est * n_sym);

figure; plot(real(Zcfo), imag(Zcfo), 'o'); axis equal; grid on
title('Part 4: After CFO correction (still unequalized)');
xlabel('Re'); ylabel('Im');

% 4.10) LS estimate of M-tap discrete equivalent channel USING CFO-corrected data
y_ls = Zcfo(M:N_tr);
h_d_ls = (train_matr' * train_matr) \ (train_matr' * y_ls);

resid = y_ls - train_matr*h_d_ls;
fprintf('[Part 4.10] LS residual mean |resid|^2 = %.3e\n', mean(abs(resid).^2));

% 4.11) ZF equalizer (L = 5M, delta = M)
h_d = h_d_ls(:);
Mch = length(h_d);
L = 5*Mch;
delta = Mch; % 0-based

H = toeplitz([h_d; zeros(L-1,1)], [h_d(1); zeros(L-1,1)]);
e_delta = zeros(Mch+L-1,1);  e_delta(delta+1) = 1;

f_zf = (H' * H) \ (H' * e_delta);
g_overall = conv(h_d, f_zf);
[pk, ipk] = max(abs(g_overall));

fprintf('[Part 4.11] max|conv(h_d,f)|=%.3f at n=%d (0-based), target delta=%d\n', ...
        pk, ipk-1, delta);

figure; stem(0:length(g_overall)-1, abs(g_overall)); grid on
title('Part 4: |conv(h_d, f_{ZF})|'); xlabel('n'); ylabel('|.|');

% 4.12) Filter, discard delta, keep next N 
Zeq = conv(Zcfo, f_zf);

start_idx = delta + 1;
end_idx   = start_idx + N - 1;
if length(Zeq) < end_idx
    Zeq = [Zeq; zeros(end_idx - length(Zeq),1)];
end

Z_keep = Zeq(start_idx:end_idx);

figure; plot(real(Z_keep), imag(Z_keep), 'o'); axis equal; grid on
title('Part 4: After CFO correction + ZF (discard \delta, keep N)');
xlabel('Re'); ylabel('Im');

fprintf('[Part 4] mean |Z_keep|^2 = %.4f\n', mean(abs(Z_keep).^2));

