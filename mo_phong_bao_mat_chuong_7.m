clc; clear; close all;

%% ===================== ENV (Al-Hourani) =====================
env = {'Suburban','Urban','Dense Urban'};
a_env = [4.88  9.61  12.08];
b_env = [0.29  0.16  0.11];     % (Al-Hourani)
EENV = numel(env);

%% ===================== SIM PARAMETERS =====================
rng(1);

Ps_dB = -40:2:50;              % dB ratio (x-axis you want)
Ps    = 10.^(Ps_dB/10);

Rs = 0.2;                       % secrecy target rate
N0 = 1e-9;                      % noise power

% Two-loop MC
Npos = 400;                     % outer: number of Eve positions
Nfad = 400;                     % inner: fading samples per position
N    = Npos*Nfad;

%% ===================== GEOMETRY =====================
H = 100;

% S = (0,0,H), Bob at horizontal distance 50
xS = 0; yS = 0; zS = H;
xB = 30; yB = 0; zB = 1.5;

% Eve region center and radius
x_cE = 120; y_cE = 0; zE = 1.5;
R_E  = 70;

% Fixed distances for Bob
r_SB_2D = sqrt((xB-xS)^2 + (yB-yS)^2);
d_SB_3D = sqrt(r_SB_2D^2 + (zS - zB)^2);
theta_SB = atan2(zS - zB, r_SB_2D) * 180/pi;

% FIX Eve at center
r_SE0_2D = sqrt((x_cE-xS)^2 + (y_cE-yS)^2);
d_SE0_3D = sqrt(r_SE0_2D^2 + (zS - zE)^2);
theta_SE0 = atan2(zS - zE, r_SE0_2D) * 180/pi;

%% ===================== PATHLOSS: FSPL + excess loss =====================
fc = 2e9;
c0 = 3e8;
eta_LoS_dB  = 1;
eta_NLoS_dB = 20;

FSPL_const_dB = 20*log10(4*pi*fc/c0);
Lfspl_SB_dB   = FSPL_const_dB + 20*log10(d_SB_3D);
Lfspl_SE0_dB  = FSPL_const_dB + 20*log10(d_SE0_3D);

%% ===================== FADING: Hybrid =====================
Kenv_dB = [10 6 2];
Kenv    = 10.^(Kenv_dB/10);

%% ===================== OUTPUT =====================
% type = 1: FIX (E = center)
% type = 2: RANDOM (E in disk)
T = 2;
Cs_avg = zeros(EENV, T, numel(Ps));
POS    = zeros(EENV, T, numel(Ps)); % SOP = P(Cs < Rs)

%% ===================== MAIN LOOP over environments =====================
for e = 1:EENV

    a = a_env(e);
    b = b_env(e);
    K = Kenv(e);

    % PLoS for Bob (fixed geometry)
    P_LOS_SB = 1/(1 + a*exp(-b*(theta_SB - a)));

    % PLoS for FIX Eve at center (fixed geometry)
    P_LOS_SE0 = 1/(1 + a*exp(-b*(theta_SE0 - a)));

    fprintf('ENV=%s | a=%.2f b=%.2f | PLoS_SB=%.4f | PLoS_SE0=%.4f | K=%.1f (%.1f dB)\n',...
        env{e}, a, b, P_LOS_SB, P_LOS_SE0, K, Kenv_dB(e));

    %% Sweep Ps
    for k = 1:numel(Ps)

        Pt = Ps(k);

        %% ===================== (1) FIX: Eve fixed at center =====================
        sumCs_fix = 0;
        out_fix   = 0;

        for p = 1:Npos

            % Per-sample LoS/NLoS + fading (inner)
            isLoS_B = (rand(1,Nfad) <= P_LOS_SB);
            isLoS_E = (rand(1,Nfad) <= P_LOS_SE0);

            % Pathloss dB per sample
            L_SB_dB = Lfspl_SB_dB + eta_NLoS_dB*ones(1,Nfad);
            L_SB_dB(isLoS_B) = Lfspl_SB_dB + eta_LoS_dB;

            L_SE_dB = Lfspl_SE0_dB + eta_NLoS_dB*ones(1,Nfad);
            L_SE_dB(isLoS_E) = Lfspl_SE0_dB + eta_LoS_dB;

            L_SB_lin = 10.^(L_SB_dB/10);
            L_SE_lin = 10.^(L_SE_dB/10);

            % Fading: Rayleigh baseline
            gB = (randn(1,Nfad)+1j*randn(1,Nfad))/sqrt(2);
            gE = (randn(1,Nfad)+1j*randn(1,Nfad))/sqrt(2);
            hB = gB; hE = gE;

            % LoS -> Rician(K) with random phase
            idxB = find(isLoS_B);
            if ~isempty(idxB)
                phiB = 2*pi*rand(1,numel(idxB));
                hLOS_B = exp(1j*phiB);
                gBL = (randn(1,numel(idxB))+1j*randn(1,numel(idxB)))/sqrt(2);
                hB(idxB) = sqrt(K/(K+1))*hLOS_B + sqrt(1/(K+1))*gBL;
            end

            idxE = find(isLoS_E);
            if ~isempty(idxE)
                phiE = 2*pi*rand(1,numel(idxE));
                hLOS_E = exp(1j*phiE);
                gEL = (randn(1,numel(idxE))+1j*randn(1,numel(idxE)))/sqrt(2);
                hE(idxE) = sqrt(K/(K+1))*hLOS_E + sqrt(1/(K+1))*gEL;
            end

            gamma_B = (Pt./L_SB_lin) .* (abs(hB).^2) ./ N0;
            gamma_E = (Pt./L_SE_lin) .* (abs(hE).^2) ./ N0;

            Cs = max(log2(1+gamma_B) - log2(1+gamma_E), 0);

            sumCs_fix = sumCs_fix + sum(Cs);
            out_fix   = out_fix + sum(Cs < Rs);
        end

        Cs_avg(e,1,k) = sumCs_fix / N;
        POS(e,1,k)    = out_fix / N;

        %% ===================== (2) RANDOM: Eve random in disk (two-loop MC) =====================
        sumCs = 0;
        outCount = 0;

        for p = 1:Npos

            % random Eve position uniform in disk
            rr  = R_E*sqrt(rand);
            phi = 2*pi*rand;
            xE  = x_cE + rr*cos(phi);
            yE  = y_cE + rr*sin(phi);

            % Eve geometry
            r_SE_2D = sqrt((xE-xS)^2 + (yE-yS)^2);
            d_SE_3D = sqrt(r_SE_2D^2 + (zS - zE)^2);
            theta_SE = atan2(zS - zE, r_SE_2D) * 180/pi;

            % PLoS for this Eve position
            P_LOS_SE = 1/(1 + a*exp(-b*(theta_SE - a)));

            % FSPL for this Eve position
            Lfspl_SE_dB = FSPL_const_dB + 20*log10(d_SE_3D);

            % Inner: per-sample LoS/NLoS + fading
            isLoS_B = (rand(1,Nfad) <= P_LOS_SB);
            isLoS_E = (rand(1,Nfad) <= P_LOS_SE);

            L_SB_dB = Lfspl_SB_dB + eta_NLoS_dB*ones(1,Nfad);
            L_SB_dB(isLoS_B) = Lfspl_SB_dB + eta_LoS_dB;

            L_SE_dB = Lfspl_SE_dB + eta_NLoS_dB*ones(1,Nfad);
            L_SE_dB(isLoS_E) = Lfspl_SE_dB + eta_LoS_dB;

            L_SB_lin = 10.^(L_SB_dB/10);
            L_SE_lin = 10.^(L_SE_dB/10);

            gB = (randn(1,Nfad)+1j*randn(1,Nfad))/sqrt(2);
            gE = (randn(1,Nfad)+1j*randn(1,Nfad))/sqrt(2);
            hB = gB;  hE = gE;

            idxB = find(isLoS_B);
            if ~isempty(idxB)
                phiB = 2*pi*rand(1,numel(idxB));
                hLOS_B = exp(1j*phiB);
                gBL = (randn(1,numel(idxB))+1j*randn(1,numel(idxB)))/sqrt(2);
                hB(idxB) = sqrt(K/(K+1))*hLOS_B + sqrt(1/(K+1))*gBL;
            end

            idxE = find(isLoS_E);
            if ~isempty(idxE)
                phiE = 2*pi*rand(1,numel(idxE));
                hLOS_E = exp(1j*phiE);
                gEL = (randn(1,numel(idxE))+1j*randn(1,numel(idxE)))/sqrt(2);
                hE(idxE) = sqrt(K/(K+1))*hLOS_E + sqrt(1/(K+1))*gEL;
            end

            gamma_B = (Pt./L_SB_lin) .* (abs(hB).^2) ./ N0;
            gamma_E = (Pt./L_SE_lin) .* (abs(hE).^2) ./ N0;

            Cs = max(log2(1+gamma_B) - log2(1+gamma_E), 0);

            sumCs = sumCs + sum(Cs);
            outCount = outCount + sum(Cs < Rs);
        end

        Cs_avg(e,2,k) = sumCs / N;
        POS(e,2,k)    = outCount / N;

    end % end Ps sweep
end % end env

%% ===================== PLOTS (x-axis = Ps_dB, TOP-DOWN layout) =====================
type_names = {'FIX: E=c_E', 'RANDOM: E in disk'};
lineStyles = {'-','--'};         % FIX solid, RANDOM dashed

figure('Color','w','Name','STAGE 2 | FIX vs RANDOM | Multi-Environment | x=Ps(dB)');

mainTitle = sprintf(['STAGE 2 | FIX vs RANDOM | Prob LoS/NLoS (Al-Hourani) | ',...
    'FSPL+excess (etaL=%gdB, etaN=%gdB) | Npos=%d, Nfad=%d | Rs=%.2f | x=Ps(dB)'], ...
    eta_LoS_dB, eta_NLoS_dB, Npos, Nfad, Rs);
annotation(gcf,'textbox',[0 0.95 1 0.05],'String',mainTitle, ...
    'EdgeColor','none','HorizontalAlignment','center','FontWeight','bold');

% legend: each env has 2 curves
leg = cell(1, 2*EENV);
ii = 1;
for e = 1:EENV
    leg{ii} = sprintf('%s | %s', env{e}, type_names{1}); ii = ii + 1;
    leg{ii} = sprintf('%s | %s', env{e}, type_names{2}); ii = ii + 1;
end

% ---------- TOP: Avg secrecy rate ----------
subplot(2,1,1); hold on; grid on;
for e = 1:EENV
    y_fix = squeeze(Cs_avg(e,1,:))';
    y_rnd = squeeze(Cs_avg(e,2,:))';

    plot(Ps_dB, y_fix, 'LineWidth', 2, 'LineStyle', lineStyles{1});
    plot(Ps_dB, y_rnd, 'LineWidth', 2, 'LineStyle', lineStyles{2});
end
xlabel('P_s (dB) (your sweep)');
ylabel('Avg Secrecy Rate  E[C_s]');
title('Average Secrecy Rate');
legend(leg,'Location','northwest');

% ---------- BOTTOM: SOP ----------
subplot(2,1,2); hold on; grid on;
for e = 1:EENV
    y_fix = squeeze(POS(e,1,:))';
    y_rnd = squeeze(POS(e,2,:))';

    semilogy(Ps_dB, y_fix, 'LineWidth', 2, 'LineStyle', lineStyles{1});
    semilogy(Ps_dB, y_rnd, 'LineWidth', 2, 'LineStyle', lineStyles{2});
end
xlabel('P_s (dB) (your sweep)');
ylabel('Secrecy Outage Probability  P(C_s<R_s)');
title(sprintf('SOP (R_s=%.2f)', Rs));
legend(leg,'Location','southwest');
ylim([1e-4 1]);

fprintf('\nDone.\n');
fprintf('- Same simulation/logic.\n');
fprintf('- Changed plotting: x-axis = Ps_dB and subplots are TOP-DOWN (2x1).\n');
