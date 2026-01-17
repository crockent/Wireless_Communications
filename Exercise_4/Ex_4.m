clear all, 
close all, 
clc;


T = 1;
over = 10;

Ts= T / over;

beta =0.4; %roll-off factor
B = 3; %half symbol length
Ntr = 50; 
N = 200;
%% 1
%1.1

data_bits = randi([0 1], 2*(N-Ntr), 1);
training_bits = randi([0 1], 2*Ntr, 1);

g(:,1) =  srrc_pulse(T,Ts,B,beta);

M =10;

data_s = bits_to_4qam(data_bits);
trainig_s = bits_to_4qam(training_bits);

packet=[trainig_s;data_s];

packet = upsample(packet,over);

A = conv(packet,g);

% 1.2

c = 1; %c = (randn + 1i*randn)/sqrt(2); change c if you want 1.8

Y_r = conv(A,c)*Ts; %received from channel

%1.3
Z_r = conv(Y_r,g)*Ts; %filtered 

%1.4;
analog_channel = conv(conv(c,g),g)*Ts;
x = length(analog_channel) -1;
figure;
plot(0:1:x,abs(analog_channel));
xlabel("Time");
ylabel("Amplitude")
grid on;
title("Complex Analog Channel");

%1.5
[Ed,d_opt] = energy_Synchronization(Z_r, over, N,analog_channel);

figure;
plot(1:length(Ed), Ed); grid on;
title('E_d '); 
xlabel('d (samples)'); ylabel('E_d');

% Σύγκριση με composite channel
figure; hold on;
plot(real(analog_channel)); plot(imag(analog_channel));
sample_idx = d_opt : over : d_opt + (Ntr-1)*over;
sample_idx = sample_idx(sample_idx <= length(analog_channel));
stem(sample_idx, real(analog_channel(sample_idx)), '--b');
stem(sample_idx, imag(analog_channel(sample_idx)), '--r');
grid on;

% 1.6

[C_d, corr_Y, d_opt_corr] = training_sync(Z_r, trainig_s, over, Ntr ,analog_channel);

s_en = 2; %symbol energy

h_hat = C_d / (Ntr*(s_en)^2);

d_total = length(analog_channel) - B * over;

% Plot 2: Normalized |C_d| vs |analog_channel| (composite)
Lplot = min(d_total, length(analog_channel));
figure; hold on; grid on;
plot((corr_Y(1:Lplot)/max(corr_Y(1:Lplot))), 'LineWidth', 1.2);
plot(abs(analog_channel(1:Lplot))/max(abs(analog_channel(1:Lplot))), '--', 'LineWidth', 1.2);
title('Training metric vs |h_{composite}| (norm)');
xlabel('d (samples)'); ylabel('normalized');

%1.7





%% 2
clear all, 
close all, 
clc;

T = 1;
over = 10;

Ts= T / over;

beta =0.4; %roll-off factor
B = 3; %half symbol length
Ntr = 50; 
N = 200;

%2.1

data_bits = randi([0 1], 2*(N-Ntr), 1);
training_bits = randi([0 1], 2*Ntr, 1);

g(:,1) =  srrc_pulse(T,Ts,B,beta);

M =10;

data_s = bits_to_4qam(data_bits);
trainig_s = bits_to_4qam(training_bits);

packet=[trainig_s;data_s];

packet = upsample(packet,over);

A = conv(packet,g);

c = 1; %c = (randn + 1i*randn)/sqrt(2); change c if you want 1.8 step

DeltaF = 10^(-3);
phi = 0;
Y_r = conv(A,c) * Ts;
t = (0:length(Y_r)-1).';
Y_r_CFO = conv(A,c).*exp(1j * (2*pi*DeltaF *t )); %received from channel

%2.2
Z_r = conv(Y_r_CFO,g)*Ts; %filtered 

analog_channel = conv(conv(c,g),g)*Ts;
x = length(analog_channel) -1;
figure;
plot(0:1:x,abs(analog_channel));
xlabel("Time");
ylabel("Amplitude")
grid on;
title("Complex Analog Channel");


[Ed,d_opt] = energy_Synchronization(Z_r, over, N,analog_channel);

figure;
plot(1:length(Ed), Ed); grid on;
title('E_d '); 
xlabel('d (samples)'); ylabel('E_d');

% Σύγκριση με composite channel
figure; hold on;
plot(real(analog_channel)); plot(imag(analog_channel));
sample_idx = d_opt : over : d_opt + (Ntr-1)*over;
sample_idx = sample_idx(sample_idx <= length(analog_channel));
stem(sample_idx, real(analog_channel(sample_idx)), '--b');
stem(sample_idx, imag(analog_channel(sample_idx)), '--r');
grid on;

[C_d, corr_Y, d_opt_corr] = training_sync(Z_r, trainig_s, over, Ntr ,analog_channel);

s_en = 2; %symbol energy

h_hat = C_d / (Ntr*(s_en)^2);

d_total = length(analog_channel) - B * over;

% Plot 2: Normalized |C_d| vs |analog_channel| (composite)
Lplot = min(d_total, length(analog_channel));
figure; hold on; grid on;
plot((corr_Y(1:Lplot)/max(corr_Y(1:Lplot))), 'LineWidth', 1.2);
plot(abs(analog_channel(1:Lplot))/max(abs(analog_channel(1:Lplot))), '--', 'LineWidth', 1.2);
title('Training metric vs |h_{composite}| (norm)');
xlabel('d (samples)'); ylabel('normalized');

%1.7
