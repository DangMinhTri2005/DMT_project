clc; clear; close all;

%% ===================== ONLY URBAN (Al-Hourani) =====================
env_name = 'Urban';
a = 9.61;
b = 0.16;

%% ===================== SIM PARAMETERS =====================
rng(1);

Ps_dB = -40:2:100;              % dB ratio  (X-AXIS)
Ps    = 10.^(Ps_dB/10);         % linear ratio (consistent with your code)

Rs = 0.2;                       % secrecy target rate

% IMPORTANT: reasonable wireless scale for FSPL
N0 = 1e-9;                      % noise power (W)

% Two-loop MC
Npos = 400;                     % outer: number of Eve positions
Nfad = 400;                     % inner: fading samples per position
N    = Npos*Nfad;

%% ===================== GEOMETRY =====================
H = 100;

% S = (0,0,H), Bob
xS = 0; yS = 0; zS = H;
xB = 30; yB = 0; zB = 1.5;

% Eve region center and radius
x_cE = -120; y_cE = 0; zE = 1.5;
R_E  = 20;

% Friendly jammer placement
xJ = x_cE + 5;
yJ = y_cE + 5;
zJ = H-30;

% Fixed distances for Bob
r_SB_2D  = sqrt((xB-xS)^2 + (yB-yS)^2);
d_SB_3D  = sqrt(r_SB_2D^2 + (zS - zB)^2);
theta_SB = atan2(zS - zB, r_SB_2D) * 180/pi;

% FIX Eve at center
r_SE0_2D  = sqrt((x_cE-xS)^2 + (y_cE-yS)^2);
d_SE0_3D  = sqrt(r_SE0_2D^2 + (zS - zE)^2);
theta_SE0 = atan2(zS - zE, r_SE0_2D) * 180/pi;

% Jammer -> Bob (fixed)
r_JB_2D = sqrt((xJ-xB)^2 + (yJ-yB)^2);
d_JB_3D = sqrt(r_JB_2D^2 + (zJ - zB)^2);
theta_JB = atan2(zJ - zB, r_JB_2D) * 180/pi;

% Jammer -> Eve center (fixed for FIX)
r_JE0_2D = sqrt((xJ-x_cE)^2 + (yJ-y_cE)^2);
d_JE0_3D = sqrt(r_JE0_2D^2 + (zJ - zE)^2);
theta_JE0 = atan2(zJ - zE, r_JE0_2D) * 180/pi;

%% ===================== PATHLOSS: FSPL + excess loss =====================
fc = 2e9;
c0 = 3e8;
eta_LoS_dB  = 1;
eta_NLoS_dB = 20;

FSPL_const_dB = 20*log10(4*pi*fc/c0);

Lfspl_SB_dB  = FSPL_const_dB + 20*log10(d_SB_3D);
Lfspl_SE0_dB = FSPL_const_dB + 20*log10(d_SE0_3D);
Lfspl_JB_dB  = FSPL_const_dB + 20*log10(d_JB_3D);
Lfspl_JE0_dB = FSPL_const_dB + 20*log10(d_JE0_3D);

%% ===================== FADING: Hybrid (LoS Rician, NLoS Rayleigh) =====================
K_dB = 6;
K    = 10^(K_dB/10);

%% ===================== LoS probabilities for fixed links (S-links + J-links) =====================
P_LOS_SB  = 1/(1 + a*exp(-b*(theta_SB  - a)));
P_LOS_SE0 = 1/(1 + a*exp(-b*(theta_SE0 - a)));

% NEW (Cách 2): PLoS for jammer links (fixed geometry)
P_LOS_JB  = 1/(1 + a*exp(-b*(theta_JB  - a)));
P_LOS_JE0 = 1/(1 + a*exp(-b*(theta_JE0 - a)));

fprintf('ENV=%s | a=%.2f b=%.2f | K=%g dB\n', env_name, a, b, K_dB);
fprintf('P_LOS_SB=%.4f | P_LOS_SE0=%.4f\n', P_LOS_SB, P_LOS_SE0);
fprintf('P_LOS_JB=%.4f | P_LOS_JE0=%.4f\n', P_LOS_JB, P_LOS_JE0);
fprintf('J=(%.1f,%.1f,%.1f) | d_JB=%.2f | d_JE0=%.2f\n', xJ,yJ,zJ,d_JB_3D,d_JE0_3D);

%% ===================== Jammer power settings =====================
PJ_dBm_list = [-Inf, 40];       % NoJam and Jam=40 dBm (10 W)
PJ_list = zeros(size(PJ_dBm_list));
for i = 1:numel(PJ_dBm_list)
    if isinf(PJ_dBm_list(i))
        PJ_list(i) = 0;
    else
        PJ_list(i) = 10.^((PJ_dBm_list(i) - 30)/10);
    end
end
nPJ = numel(PJ_list);

%% ===================== OUTPUT =====================
% type = 1: FIX (E center)
% type = 2: RANDOM (E in disk)
T = 2;

Cs_avg = zeros(T, nPJ, numel(Ps));
POS    = zeros(T, nPJ, numel(Ps));

%% ===================== MAIN LOOP over Ps =====================
for k = 1:numel(Ps)

    Pt = Ps(k);

    %% ===================== (1) FIX: Eve fixed at center =====================
    % CRN: share SAME random samples for NoJam/Jam to reduce variance
    sumCs_fix = zeros(1,nPJ);
    out_fix   = zeros(1,nPJ);

    for p = 1:Npos

        % Shared randomness for S-links
        isLoS_B  = (rand(1,Nfad) <= P_LOS_SB);
        isLoS_E0 = (rand(1,Nfad) <= P_LOS_SE0);

        % NEW (Cách 2): Shared randomness for J-links (fixed geometry)
        isLoS_JB  = (rand(1,Nfad) <= P_LOS_JB);
        isLoS_JE0 = (rand(1,Nfad) <= P_LOS_JE0);

        % Pathloss S->B (per-sample)
        L_SB_dB = Lfspl_SB_dB + eta_NLoS_dB*ones(1,Nfad);
        L_SB_dB(isLoS_B) = Lfspl_SB_dB + eta_LoS_dB;

        % Pathloss S->E0 (per-sample)
        L_SE0_dB = Lfspl_SE0_dB + eta_NLoS_dB*ones(1,Nfad);
        L_SE0_dB(isLoS_E0) = Lfspl_SE0_dB + eta_LoS_dB;

        L_SB_lin  = 10.^(L_SB_dB/10);
        L_SE0_lin = 10.^(L_SE0_dB/10);

        % NEW (Cách 2): Pathloss J->B, J->E0 (per-sample, Al-Hourani LoS/NLoS)
        L_JB_dB  = Lfspl_JB_dB  + eta_NLoS_dB*ones(1,Nfad);
        L_JB_dB(isLoS_JB) = Lfspl_JB_dB + eta_LoS_dB;

        L_JE0_dB = Lfspl_JE0_dB + eta_NLoS_dB*ones(1,Nfad);
        L_JE0_dB(isLoS_JE0) = Lfspl_JE0_dB + eta_LoS_dB;

        L_JB_lin  = 10.^(L_JB_dB/10);
        L_JE0_lin = 10.^(L_JE0_dB/10);

        % Fading baseline for S-links
        gB  = (randn(1,Nfad)+1j*randn(1,Nfad))/sqrt(2);
        gE0 = (randn(1,Nfad)+1j*randn(1,Nfad))/sqrt(2);
        hB  = gB;
        hE0 = gE0;

        % LoS -> Rician(K) for S->B
        idxB = find(isLoS_B);
        if ~isempty(idxB)
            phiB = 2*pi*rand(1,numel(idxB));
            hLOS_B = exp(1j*phiB);
            gBL = (randn(1,numel(idxB))+1j*randn(1,numel(idxB)))/sqrt(2);
            hB(idxB) = sqrt(K/(K+1))*hLOS_B + sqrt(1/(K+1))*gBL;
        end

        % LoS -> Rician(K) for S->E0
        idxE0 = find(isLoS_E0);
        if ~isempty(idxE0)
            phiE = 2*pi*rand(1,numel(idxE0));
            hLOS_E = exp(1j*phiE);
            gEL = (randn(1,numel(idxE0))+1j*randn(1,numel(idxE0)))/sqrt(2);
            hE0(idxE0) = sqrt(K/(K+1))*hLOS_E + sqrt(1/(K+1))*gEL;
        end

        % Jammer fading (Rayleigh)
        hJB  = (randn(1,Nfad)+1j*randn(1,Nfad))/sqrt(2);
        hJE0 = (randn(1,Nfad)+1j*randn(1,Nfad))/sqrt(2);

        % Evaluate both NoJam/Jam on same samples
        for jp = 1:nPJ
            PJ = PJ_list(jp);

            SINR_B = (Pt./L_SB_lin)  .* (abs(hB).^2)  ./ (N0 + (PJ./L_JB_lin) .* (abs(hJB).^2));
            SINR_E = (Pt./L_SE0_lin) .* (abs(hE0).^2) ./ (N0 + (PJ./L_JE0_lin).* (abs(hJE0).^2));

            Cs = max(log2(1+SINR_B) - log2(1+SINR_E), 0);

            sumCs_fix(jp) = sumCs_fix(jp) + sum(Cs);
            out_fix(jp)   = out_fix(jp)   + sum(Cs < Rs);
        end
    end

    for jp = 1:nPJ
        Cs_avg(1,jp,k) = sumCs_fix(jp)/N;
        POS(1,jp,k)    = out_fix(jp)/N;
    end

    %% ===================== (2) RANDOM: Eve random in disk (two-loop MC) =====================
    sumCs_rnd = zeros(1,nPJ);
    out_rnd   = zeros(1,nPJ);

    for p = 1:Npos

        % random Eve position uniform in disk
        rr  = R_E*sqrt(rand);
        ph  = 2*pi*rand;
        xE  = x_cE + rr*cos(ph);
        yE  = y_cE + rr*sin(ph);

        % Eve geometry
        r_SE_2D  = sqrt((xE-xS)^2 + (yE-yS)^2);
        d_SE_3D  = sqrt(r_SE_2D^2 + (zS - zE)^2);
        theta_SE = atan2(zS - zE, r_SE_2D) * 180/pi;

        % PLoS for this Eve position (S->E)
        P_LOS_SE = 1/(1 + a*exp(-b*(theta_SE - a)));

        % FSPL for this Eve position (S->E)
        Lfspl_SE_dB = FSPL_const_dB + 20*log10(d_SE_3D);

        % Jammer->E distance for this Eve position
        r_JE_2D = sqrt((xJ-xE)^2 + (yJ-yE)^2);
        d_JE_3D = sqrt(r_JE_2D^2 + (zJ - zE)^2);
        theta_JE = atan2(zJ - zE, r_JE_2D) * 180/pi;

        % NEW (Cách 2): PLoS for jammer->E at this position
        P_LOS_JE = 1/(1 + a*exp(-b*(theta_JE - a)));

        % FSPL for jammer->E at this position
        Lfspl_JE_dB = FSPL_const_dB + 20*log10(d_JE_3D);

        % Shared inner randomness (S-links + J-links)
        isLoS_B  = (rand(1,Nfad) <= P_LOS_SB);
        isLoS_E  = (rand(1,Nfad) <= P_LOS_SE);

        % NEW (Cách 2): jammer link states
        isLoS_JB = (rand(1,Nfad) <= P_LOS_JB);
        isLoS_JE = (rand(1,Nfad) <= P_LOS_JE);

        % S->B loss
        L_SB_dB = Lfspl_SB_dB + eta_NLoS_dB*ones(1,Nfad);
        L_SB_dB(isLoS_B) = Lfspl_SB_dB + eta_LoS_dB;

        % S->E loss
        L_SE_dB = Lfspl_SE_dB + eta_NLoS_dB*ones(1,Nfad);
        L_SE_dB(isLoS_E) = Lfspl_SE_dB + eta_LoS_dB;

        L_SB_lin = 10.^(L_SB_dB/10);
        L_SE_lin = 10.^(L_SE_dB/10);

        % NEW (Cách 2): J->B and J->E losses
        L_JB_dB = Lfspl_JB_dB + eta_NLoS_dB*ones(1,Nfad);
        L_JB_dB(isLoS_JB) = Lfspl_JB_dB + eta_LoS_dB;

        L_JE_dB = Lfspl_JE_dB + eta_NLoS_dB*ones(1,Nfad);
        L_JE_dB(isLoS_JE) = Lfspl_JE_dB + eta_LoS_dB;

        L_JB_lin = 10.^(L_JB_dB/10);
        L_JE_lin = 10.^(L_JE_dB/10);

        % Fading baseline
        gB = (randn(1,Nfad)+1j*randn(1,Nfad))/sqrt(2);
        gE = (randn(1,Nfad)+1j*randn(1,Nfad))/sqrt(2);
        hB = gB;
        hE = gE;

        % LoS -> Rician(K) for S->B
        idxB = find(isLoS_B);
        if ~isempty(idxB)
            phiB = 2*pi*rand(1,numel(idxB));
            hLOS_B = exp(1j*phiB);
            gBL = (randn(1,numel(idxB))+1j*randn(1,numel(idxB)))/sqrt(2);
            hB(idxB) = sqrt(K/(K+1))*hLOS_B + sqrt(1/(K+1))*gBL;
        end

        % LoS -> Rician(K) for S->E
        idxE = find(isLoS_E);
        if ~isempty(idxE)
            phiE = 2*pi*rand(1,numel(idxE));
            hLOS_E = exp(1j*phiE);
            gEL = (randn(1,numel(idxE))+1j*randn(1,numel(idxE)))/sqrt(2);
            hE(idxE) = sqrt(K/(K+1))*hLOS_E + sqrt(1/(K+1))*gEL;
        end

        % Jammer fading (Rayleigh)
        hJB = (randn(1,Nfad)+1j*randn(1,Nfad))/sqrt(2);
        hJE = (randn(1,Nfad)+1j*randn(1,Nfad))/sqrt(2);

        % Evaluate both NoJam/Jam on same samples
        for jp = 1:nPJ
            PJ = PJ_list(jp);

            SINR_B = (Pt./L_SB_lin) .* (abs(hB).^2) ./ (N0 + (PJ./L_JB_lin).* (abs(hJB).^2));
            SINR_E = (Pt./L_SE_lin) .* (abs(hE).^2) ./ (N0 + (PJ./L_JE_lin).* (abs(hJE).^2));

            Cs = max(log2(1+SINR_B) - log2(1+SINR_E), 0);

            sumCs_rnd(jp) = sumCs_rnd(jp) + sum(Cs);
            out_rnd(jp)   = out_rnd(jp)   + sum(Cs < Rs);
        end
    end

    for jp = 1:nPJ
        Cs_avg(2,jp,k) = sumCs_rnd(jp)/N;
        POS(2,jp,k)    = out_rnd(jp)/N;
    end

end % k

%% ===================== PLOTS: TOP (Cs), BOTTOM (SOP) | x = Ps(dB) =====================
type_names = {'FIX: E=c_E', 'RANDOM: E in disk'};
jam_names  = {'No Jam', sprintf('Jam PJ=%g dBm', PJ_dBm_list(2))};

leg = { ...
    sprintf('%s | %s', type_names{1}, jam_names{1}), ...
    sprintf('%s | %s', type_names{2}, jam_names{1}), ...
    sprintf('%s | %s', type_names{1}, jam_names{2}), ...
    sprintf('%s | %s', type_names{2}, jam_names{2}) ...
};

x_plot = Ps_dB;

Cs_fix_no  = squeeze(Cs_avg(1,1,:))';
Cs_rnd_no  = squeeze(Cs_avg(2,1,:))';
Cs_fix_j   = squeeze(Cs_avg(1,2,:))';
Cs_rnd_j   = squeeze(Cs_avg(2,2,:))';

SOP_fix_no = squeeze(POS(1,1,:))';
SOP_rnd_no = squeeze(POS(2,1,:))';
SOP_fix_j  = squeeze(POS(1,2,:))';
SOP_rnd_j  = squeeze(POS(2,2,:))';

figure('Color','w','Name','Urban | Friendly Jamming | FIX vs RANDOM | x=Ps(dB) | Al-Hourani for jammer links');

mainTitle = sprintf(['ENV=%s | Al-Hourani LoS/NLoS for S-links + J-links | FSPL+excess (etaL=%gdB, etaN=%gdB) | ',...
    'Hybrid fading (LoS Rician K=%gdB, NLoS Rayleigh) | J=(%.0f,%.0f,%.0f) | Rs=%.2f | Npos=%d,Nfad=%d | x=Ps(dB)'], ...
    env_name, eta_LoS_dB, eta_NLoS_dB, K_dB, xJ,yJ,zJ, Rs, Npos, Nfad);
annotation(gcf,'textbox',[0 0.95 1 0.05],'String',mainTitle, ...
    'EdgeColor','none','HorizontalAlignment','center','FontWeight','bold');

subplot(2,1,1); hold on; grid on;
plot(x_plot, Cs_fix_no, 'LineWidth', 2, 'LineStyle', '-');
plot(x_plot, Cs_rnd_no, 'LineWidth', 2, 'LineStyle', '--');
plot(x_plot, Cs_fix_j,  'LineWidth', 2, 'LineStyle', '-');
plot(x_plot, Cs_rnd_j,  'LineWidth', 2, 'LineStyle', '--');
ylabel('Avg Secrecy Rate  E[C_s]');
title('Average Secrecy Rate');
legend(leg,'Location','northwest');

subplot(2,1,2); hold on; grid on;
semilogy(x_plot, SOP_fix_no, 'LineWidth', 2, 'LineStyle', '-');
semilogy(x_plot, SOP_rnd_no, 'LineWidth', 2, 'LineStyle', '--');
semilogy(x_plot, SOP_fix_j,  'LineWidth', 2, 'LineStyle', '-');
semilogy(x_plot, SOP_rnd_j,  'LineWidth', 2, 'LineStyle', '--');
xlabel('P_s (dB) (your sweep)');
ylabel('SOP = P(C_s<R_s)');
title(sprintf('Secrecy Outage Probability (R_s=%.2f)', Rs));
legend(leg,'Location','southwest');
ylim([1e-4 1]);

fprintf('\nDone.\n');
fprintf('- Only URBAN used.\n');
fprintf('- Cách 2 applied: jammer links J->B and J->E use Al-Hourani LoS/NLoS per-sample (NOT fixed +20dB).\n');
fprintf('- CRN applied: NoJam/Jam share the SAME random samples per (Ps, position, fading), so curves are smoother.\n');
fprintf('- Plot x-axis is Ps_dB (no sorting), so no kink from re-ordering.\n');
