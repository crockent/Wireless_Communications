clear all, close all, clc;


T = 1;
over = 10;

Ts= T / over;

beta =0.4; %roll-off factor
B = 3; %half symbol length
Ntr = 50; 

%% 1
%1.1
N = 200;
data_bits = randi([0 1], 2*N-2*Ntr, 1);
training_bits = randi([0 1], 2*Ntr, 1);

g(:,1) =  srrc_pulse(T,Ts,B,beta);


data_s = bits_to_4qam(data_bits);
trainig_s = bits_to_4qam(training_bits);

packet=[trainig_s;data_s];

A = conv(packet,g);

% 1.2

c = 1;
Y_r = conv(A,c); %received from channel

%1.3
Z_r = conv(Y_r,g)*Ts; %filtered 

%1.4;
analog_channel = c * conv(g,g)*Ts;
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
[C_d, corr_Y, d_opt_corr] = training_sync(Z_r, trainig_s, over, Ntr,analog_channel);

figure;
plot(1:length(corr_Y), corr_Y); grid on;
title('Training-based metric |C_d|'); 

s_en = 2; %symbol energy

h_hat = C_d / (Ntr*(s_en)^2);

[Ed_tr,d_opt_tr] = energy_Synchronization(h_hat, over, Ntr,analog_channel);


figure;
plot(1:length(Ed_tr), Ed_tr); grid on;
title('E_d from h_hat'); 

% Plot 2: Normalized |C_d| vs |analog_channel| (composite)
Lplot = min(length(corr_Y), length(analog_channel));
figure; hold on; grid on;
plot(1:Lplot, corr_Y(1:Lplot)/max(corr_Y(1:Lplot)+eps), 'LineWidth', 1.2);
plot(1:Lplot, abs(analog_channel(1:Lplot))/max(abs(analog_channel(1:Lplot))+eps), ...
     '--', 'LineWidth', 1.2);
title('Training metric vs |h_{composite}| (norm)');
xlabel('d (samples)'); ylabel('normalized');







