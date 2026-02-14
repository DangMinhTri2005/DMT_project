% Compare Rayleigh, Rician, and Nakagami-m fading on the same plots
% Sweep Ps (dBm) -> compute instantaneous SNR gamma = (Ps/Loss)*|h|^2/N0
% Then plot Capacity / BER / Outage Probability vs Average SNR (dB)

clear; close all; clc; tic;

%% ===================== Common parameters =====================
n = 1e6;
bits = rand(1,n) > 0.5;
s = 2*bits - 1;                        % BPSK, E|s|^2=1

% Target rate for outage (bps/Hz)
R = 1;
gamma_th = 2^R - 1;

% Large-scale loss
d = 30;
alpha = 2;
Loss = d^alpha;                        % Loss > 1

% Transmit power sweep (use ONE sweep for fair comparison)
Ps_dBm = -40:2:60;
Ps = 10.^((Ps_dBm - 30)/10);           % W

% Noise power (W)
N0 = 1e-9;

%% ===================== 1) Rayleigh =====================
BER_ray = zeros(size(Ps));
OP_ray  = zeros(size(Ps));
SNRavg_ray_lin = zeros(size(Ps));
C_ray = zeros(size(Ps));

for i = 1:length(Ps)
    h = (randn(1,n) + 1j*randn(1,n)) / sqrt(2);           % E|h|^2 = 1
    noise = sqrt(N0/2) * (randn(1,n) + 1j*randn(1,n));

    y = sqrt(Ps(i)/Loss) * (h .* s) + noise;
    g_tot = sqrt(Ps(i)/Loss) * h;
    sHat = y ./ g_tot;

    bitsHat = real(sHat) > 0;
    BER_ray(i) = mean(bitsHat ~= bits);

    gamma = (Ps(i)/Loss) * (abs(h).^2) / N0;
    SNRavg_ray_lin(i) = mean(gamma);
    C_ray(i) = mean(log2(1 + gamma));
    OP_ray(i) = mean(gamma < gamma_th);
end
SNRavg_ray_dB = 10*log10(SNRavg_ray_lin);
[SNR_ray_dB, idxRay] = sort(SNRavg_ray_dB);
BER_ray = BER_ray(idxRay);
C_ray   = C_ray(idxRay);
OP_ray  = OP_ray(idxRay);

%% ===================== 2) Rician =====================
K_dB = 5;
phi = 0;
K = 10^(K_dB/10);

BER_ric = zeros(size(Ps));
OP_ric  = zeros(size(Ps));
SNRavg_ric_lin = zeros(size(Ps));
C_ric = zeros(size(Ps));

for i = 1:length(Ps)
    g = (randn(1,n) + 1j*randn(1,n)) / sqrt(2);           % CN(0,1)
    h = sqrt(K/(K+1)) * exp(1j*phi) + sqrt(1/(K+1)) * g;  % E|h|^2 = 1

    noise = sqrt(N0/2) * (randn(1,n) + 1j*randn(1,n));

    y = sqrt(Ps(i)/Loss) * (h .* s) + noise;
    g_tot = sqrt(Ps(i)/Loss) * h;
    sHat = y ./ g_tot;

    bitsHat = real(sHat) > 0;
    BER_ric(i) = mean(bitsHat ~= bits);

    gamma = (Ps(i)/Loss) * (abs(h).^2) / N0;
    SNRavg_ric_lin(i) = mean(gamma);
    C_ric(i) = mean(log2(1 + gamma));
    OP_ric(i) = mean(gamma < gamma_th);
end
SNRavg_ric_dB = 10*log10(SNRavg_ric_lin);
[SNR_ric_dB, idxRic] = sort(SNRavg_ric_dB);
BER_ric = BER_ric(idxRic);
C_ric   = C_ric(idxRic);
OP_ric  = OP_ric(idxRic);

%% ===================== 3) Nakagami-m =====================
m = 3;                                  % m>=0.5
Omega = 1;                              % E|h|^2

theta = Omega/m;                        % so E{|h|^2}=Omega

BER_nak = zeros(size(Ps));
OP_nak  = zeros(size(Ps));
SNRavg_nak_lin = zeros(size(Ps));
C_nak = zeros(size(Ps));

for i = 1:length(Ps)
    % |h|^2 ~ Gamma(m, theta), phase uniform
    gpow = gamrnd(m, theta, 1, n);
    ph = 2*pi*rand(1,n);
    h = sqrt(gpow) .* exp(1j*ph);

    noise = sqrt(N0/2) * (randn(1,n) + 1j*randn(1,n));

    y = sqrt(Ps(i)/Loss) * (h .* s) + noise;
    g_tot = sqrt(Ps(i)/Loss) * h;
    sHat = y ./ g_tot;

    bitsHat = real(sHat) > 0;
    BER_nak(i) = mean(bitsHat ~= bits);

    gamma = (Ps(i)/Loss) * (abs(h).^2) / N0;  % abs(h).^2 = gpow
    SNRavg_nak_lin(i) = mean(gamma);
    C_nak(i) = mean(log2(1 + gamma));
    OP_nak(i) = mean(gamma < gamma_th);
end
SNRavg_nak_dB = 10*log10(SNRavg_nak_lin);
[SNR_nak_dB, idxNak] = sort(SNRavg_nak_dB);
BER_nak = BER_nak(idxNak);
C_nak   = C_nak(idxNak);
OP_nak  = OP_nak(idxNak);

%% ===================== Combined plots =====================
figure;
plot(SNR_ray_dB, C_ray, 'LineWidth', 1.5); hold on; grid on;
plot(SNR_ric_dB, C_ric, 'LineWidth', 1.5);
plot(SNR_nak_dB, C_nak, 'LineWidth', 1.5);
xlabel('Average SNR (dB)'); ylabel('Average Capacity (bits/s/Hz)');
title(sprintf('Capacity vs Average SNR (R = %g bps/Hz)', R));
legend('Rayleigh', sprintf('Rician (K=%g dB)', K_dB), sprintf('Nakagami-m (m=%g)', m), 'Location', 'best');

figure;
semilogy(SNR_ray_dB, BER_ray, 'LineWidth', 1.5); hold on; grid on;
semilogy(SNR_ric_dB, BER_ric, 'LineWidth', 1.5);
semilogy(SNR_nak_dB, BER_nak, 'LineWidth', 1.5);
xlabel('Average SNR (dB)'); ylabel('BER');
title('BER vs Average SNR (BPSK)');
legend('Rayleigh', sprintf('Rician (K=%g dB)', K_dB), sprintf('Nakagami-m (m=%g)', m), 'Location', 'best');

figure;
semilogy(SNR_ray_dB, OP_ray, 'LineWidth', 1.5); hold on; grid on;
semilogy(SNR_ric_dB, OP_ric, 'LineWidth', 1.5);
semilogy(SNR_nak_dB, OP_nak, 'LineWidth', 1.5);
xlabel('Average SNR (dB)'); ylabel('Outage Probability');
title(sprintf('Outage Probability vs Average SNR (C < R, R = %g bps/Hz)', R));
legend('Rayleigh', sprintf('Rician (K=%g dB)', K_dB), sprintf('Nakagami-m (m=%g)', m), 'Location', 'best');

%% ===================== SNR vs Ps =====================
figure;
plot(Ps_dBm, 10*log10(SNRavg_ray_lin), 'LineWidth', 1.5); hold on; grid on;
plot(Ps_dBm, 10*log10(SNRavg_ric_lin), 'LineWidth', 1.5);
plot(Ps_dBm, 10*log10(SNRavg_nak_lin), 'LineWidth', 1.5);
xlabel('Transmit Power P_s (dBm)');
ylabel('Average SNR (dB)');
title('Average SNR vs Transmit Power');
legend('Rayleigh', sprintf('Rician (K=%g dB)', K_dB), sprintf('Nakagami-m (m=%g)', m), 'Location', 'best');


toc;
