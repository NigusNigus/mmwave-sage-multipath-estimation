clear; clc; close all;

%% ================== PARAMETERS ==================
ChPara = struct();

% Carrier / waveform
ChPara.fc = 28e9;
ChPara.fs = 100e6;
ChPara.Ts = 1/ChPara.fs;
ChPara.Tsig = 1e-6;

ChPara.c = 3e8;
ChPara.lambda = ChPara.c/ChPara.fc;

ChPara.SNR_dB = 20;  

% Array
ChPara.no_elements = 4;
ChPara.n = (0:ChPara.no_elements-1).';
ChPara.d = ChPara.lambda/2;

% Beamsteering Angle
ChPara.Ang_steer = -90*pi/180:10*pi/180:90*pi/180;
ChPara.thetaTX = ChPara.Ang_steer;
ChPara.thetaRX = ChPara.Ang_steer;

ChPara.Ptx = numel(ChPara.thetaTX);
ChPara.Prx = numel(ChPara.thetaRX);

% Delay grid
ChPara.delay   = 0:ChPara.Ts:ChPara.Tsig-ChPara.Ts;
ChPara.N_delay = numel(ChPara.delay);
ChPara.f = (0:ChPara.N_delay-1) * (ChPara.fs/ChPara.N_delay);

%% ================== GEOMETRY ==================
Tx = [0 0];
Rx = [0 7];

Ref_list_true = [-2.5 2.5;
                  3 3];

%% ================== BUILD TRUE MPCs ==================
[Segments_true, SegNames_true, SegMeta_true] = ...
    build_mpc_segments(Tx, Rx, Ref_list_true);

Paths_true = build_paths(Tx, Rx, Ref_list_true, Segments_true, SegMeta_true);
Npaths_true = numel(Paths_true);

ChPara.N = Npaths_true;

%% ================== SINGLE MEASUREMENT CHANNEL ==================
Measured_clean = zeros(ChPara.Ptx, ChPara.Prx, ChPara.N_delay);

for p = 1:Npaths_true

    thTX_MPC = Paths_true(p).thetaTX;
    thRX_MPC = Paths_true(p).thetaRX;

    %% --- Tx steering responses ---
    RF_vector_Tx = zeros(1, ChPara.Ptx);
    for k = 1:ChPara.Ptx
        ftx = exp(-1i*(ChPara.n) * 2*pi*ChPara.d/ChPara.lambda * ...
                  (sin(ChPara.thetaTX(k)) - sin(thTX_MPC)));
        RF_vector_Tx(k) = sum(ftx);
    end

    %% --- Rx steering responses ---
    RF_vector_Rx = zeros(1, ChPara.Prx);
    for k = 1:ChPara.Prx
        frx = exp(-1i*(ChPara.n) * 2*pi*ChPara.d/ChPara.lambda * ...
                  (sin(ChPara.thetaRX(k)) - sin(thRX_MPC)));
        RF_vector_Rx(k) = sum(frx);
    end

    %% --- Steering matrix ---
    H_steering = conj(RF_vector_Tx.') * RF_vector_Rx;   % Ptx × Prx


    %% --- Delay response ---
    delay_resp = ifft(exp(-1i*2*pi*ChPara.f * Paths_true(p).tau));
    
    %% --- Channel contribution ---
    Measured_clean = Measured_clean + Paths_true(p).alpha .* ...
        ( reshape(H_steering, ChPara.Ptx, ChPara.Prx, 1) .* ...
          reshape(delay_resp, 1, 1, []) );

    % LoS spatial component used for SNR reference (delay response is impulse-like)
    if p==1
        MPC_LoS = Paths_true(p).alpha .* ( reshape(H_steering, ChPara.Ptx, ChPara.Prx, 1)); 
    end

end


SNR_lin = 10^(ChPara.SNR_dB/10);

% Noise is added as complex AWGN with SNR defined w.r.t. the LoS peak tap power
P_sig = mean(abs(MPC_LoS(:)).^2);  % LoS average power reference
P_noise = P_sig / SNR_lin;

noise = sqrt(P_noise/2) * (randn(size(Measured_clean)) + 1i*randn(size(Measured_clean)));
Measured_MPC = Measured_clean + noise;


%% ================== RUN SAGE ==================
No_iteration = 5;
ParaEst = SAGE(Measured_MPC, No_iteration, ChPara);



%% ================== NUMERICAL RESULTS ==================

disp('--- TRUE PARAMETERS ---');
alpha_true = [Paths_true.alpha];                 % complex
alpha_true_rel = abs(alpha_true) / max(abs(alpha_true));
for p = 1:Npaths_true
    fprintf('Path %d: AoD=%.1f°, AoA=%.1f°, Delay=%.2f ns, |α|=%.3f\n', ...
        p, ...
        rad2deg(Paths_true(p).thetaTX), ...
        rad2deg(Paths_true(p).thetaRX), ...
        Paths_true(p).tau*1e9, ...
        alpha_true_rel(p));
end

disp('--- ESTIMATED PARAMETERS (relative gains) ---');

alpha_hat = abs(ParaEst.mx_MPCE(1,:));
alpha_rel = alpha_hat / max(alpha_hat);   % normalize to strongest MPC

for l = 1:ChPara.N
    fprintf('MPC %d: AoD=%.1f°, AoA=%.1f°, Delay=%.2f ns, |α|_rel=%.3f\n', ...
        l, ...
        rad2deg(ParaEst.thetaTX_MPCE(1,l)), ...
        rad2deg(ParaEst.thetaRX_MPCE(1,l)), ...
        ParaEst.DelayE(1,l)*ChPara.Ts*1e9, ...
        alpha_rel(l));
end



%% ================== FINAL PLOTS ==================

fsz = 12;
lw  = 2.0;
ms  = 15;

% Consistent colors per path
nP = Npaths_true;
C = lines(max(nP,3));

%% ==================================================
%% FIGURE 1: Geometry (Tx, Rx, Reflectors, Paths)
%% ==================================================
figure; hold on; grid on; axis equal;
set(gca,'FontSize',fsz,'LineWidth',1.3,'FontWeight','bold');

% Nodes
plot(Tx(1),Tx(2),'rs','MarkerFaceColor','r','MarkerSize',9);
plot(Rx(1),Rx(2),'bs','MarkerFaceColor','b','MarkerSize',9);

% Reflectors
for r = 1:size(Ref_list_true,1)
    plot(Ref_list_true(r,1),Ref_list_true(r,2),'ko', ...
        'MarkerFaceColor',[.5 .5 .5],'MarkerSize',7);
end

% Paths
hPath = gobjects(nP,1);
for p = 1:nP
    col = C(p,:);
    if strcmpi(Paths_true(p).type,'LOS')
        hPath(p) = plot([Tx(1) Rx(1)], [Tx(2) Rx(2)], '-', ...
            'Color',col,'LineWidth',lw);
    else
        Ref = Ref_list_true(Paths_true(p).refIdx,:);
        hPath(p) = plot([Tx(1) Ref(1) Rx(1)], [Tx(2) Ref(2) Rx(2)], '--', ...
            'Color',col,'LineWidth',lw);
    end
end

xlabel('x (m)','FontWeight','bold');
ylabel('y (m)','FontWeight','bold');
title('Geometry with LoS and Reflected Paths','FontWeight','bold');

hTx  = plot(nan,nan,'rs','MarkerFaceColor','r','MarkerSize',9);
hRx  = plot(nan,nan,'bs','MarkerFaceColor','b','MarkerSize',9);
hRef = plot(nan,nan,'ko','MarkerFaceColor',[.5 .5 .5],'MarkerSize',7);

legend([hTx hRx hRef hPath(:).'], ...
       ['Tx','Rx','Reflector', ...
        arrayfun(@(k)sprintf('Path %d',k),1:nP,'UniformOutput',false)], ...
       'Location','northeastoutside');

%% ==================================================
%% FIGURE 2: Angle-Domain Channel Power (dB)
%% ==================================================
Pwr_ang = squeeze(mean(abs(Measured_MPC).^2,3));    % avg over delay
Pwr_dB  = 10*log10(Pwr_ang + 1e-12);

figure;
imagesc(ChPara.thetaRX*180/pi, ChPara.thetaTX*180/pi, Pwr_dB);
axis xy tight; colorbar;
set(gca,'FontSize',fsz,'LineWidth',1.3,'FontWeight','bold');
xlabel('AoA (deg)','FontWeight','bold');
ylabel('AoD (deg)','FontWeight','bold');
title('Angle-Domain Channel Power (dB)','FontWeight','bold');
hold on;

% True & Estimated markers (same color per path)
for p = 1:nP
    col = C(p,:);
    plot(rad2deg(Paths_true(p).thetaRX), ...
         rad2deg(Paths_true(p).thetaTX), ...
         'o','MarkerSize',ms,'LineWidth',2, ...
         'MarkerEdgeColor',col);
    if p <= ChPara.N
        plot(rad2deg(ParaEst.thetaRX_MPCE(1,p)), ...
             rad2deg(ParaEst.thetaTX_MPCE(1,p)), ...
             'x','MarkerSize',ms,'LineWidth',2,'Color',col);
    end
end

hTrue = plot(nan,nan,'ko','MarkerSize',ms,'LineWidth',2);
hEst  = plot(nan,nan,'kx','MarkerSize',ms,'LineWidth',2);
legend([hTrue hEst],{'True MPC','Estimated MPC'}, ...
       'Location','northeastoutside');

%% ==================================================
%% FIGURE 3: Path Gains (Alpha) — Relative, dB
%% ==================================================
alpha_true = abs([Paths_true.alpha]);
alpha_true_dB = 20*log10(alpha_true / max(alpha_true));

alpha_est = abs(ParaEst.mx_MPCE(1,:));
alpha_est_dB = 20*log10(alpha_est / max(alpha_est));

figure; hold on; grid on;
set(gca,'FontSize',fsz,'LineWidth',1.3,'FontWeight','bold');

for p = 1:nP
    col = C(p,:);
    stem(p, alpha_true_dB(p), 'LineWidth',lw, ...
        'Marker','o','MarkerSize',7,'Color',col);
    if p <= numel(alpha_est_dB)
        stem(p, alpha_est_dB(p), 'LineWidth',lw, ...
            'Marker','x','MarkerSize',8,'Color',col);
    end
end

xlabel('MPC index','FontWeight','bold');
ylabel('|α| (dB, relative)','FontWeight','bold');
title('Path Gain Comparison (Relative dB)','FontWeight','bold');

hTr = plot(nan,nan,'o','LineWidth',lw,'Color','k');
hEs = plot(nan,nan,'x','LineWidth',lw,'Color','k');
legend([hTr hEs],{'True','Estimated'},'Location','northeast');
