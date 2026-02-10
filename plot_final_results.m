function plot_final_results(Tx, Rx, Ref_list_true, Paths_true, ...
                            Measured_MPC, ParaEst, ChPara, Npaths_true)
%PLOT_FINAL_RESULTS  Generates final figures 
%
% Inputs:
%   Tx, Rx               : 1x2 coordinates
%   Ref_list_true        : Nr x 2 reflector positions
%   Paths_true           : struct array with fields (type, refIdx, thetaTX, thetaRX, alpha)
%   Measured_MPC         : Ptx x Prx x Ndelay measurement tensor
%   ParaEst              : struct with fields thetaTX_MPCE, thetaRX_MPCE, mx_MPCE
%   ChPara               : struct with fields thetaTX, thetaRX, N
%   Npaths_true          : number of true paths (for consistent coloring)

fsz = 12;
lw  = 2.0;
ms  = 15;

% Consistent colors per path
nP = Npaths_true;
C = lines(max(nP,3));


%% FIGURE 1: Geometry (Tx, Rx, Reflectors, Paths)

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



%% FIGURE 2: Beamsteering Channel Power (dB) - smooth imagesc

Pwr_ang = squeeze(mean(abs(Measured_MPC).^2,3));    % avg over delay
Pwr_dB  = 10*log10(Pwr_ang + 1e-12);

% Steering-angle axes (deg)
x = rad2deg(ChPara.thetaRX);   % Rx steering angle in degrees (x-axis)
y = rad2deg(ChPara.thetaTX);   % Tx steering angle in degrees (y-axis)

% Build original grid
[X, Y] = meshgrid(x, y);

% Finer grid for smooth visualization
xq = linspace(min(x), max(x), 2*numel(x));
yq = linspace(min(y), max(y), 2*numel(y));
[Xq, Yq] = meshgrid(xq, yq);

% Interpolate (smooth)
Pwr_dB_smooth = interp2(X, Y, Pwr_dB, Xq, Yq, 'spline');

figure;
imagesc(xq, yq, Pwr_dB_smooth);
axis xy tight; colorbar;

set(gca,'FontSize',fsz,'LineWidth',1.3,'FontWeight','bold');
xlabel('Steering Angle Rx (Deg)','FontWeight','bold');
ylabel('Steering Angle Tx (Deg)','FontWeight','bold');
title('Beamsteering-Domain Channel Power (dB)','FontWeight','bold');
hold on;

% True & Estimated markers (same color per path)
for p = 1:nP
    col = C(p,:);

    % True MPC marker
    plot(rad2deg(Paths_true(p).thetaRX), ...
         rad2deg(Paths_true(p).thetaTX), ...
         'o','MarkerSize',ms,'LineWidth',2, ...
         'MarkerEdgeColor',col);

    % Estimated MPC marker
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


%% FIGURE 3: Path Gains (Alpha) — Relative, dB

alpha_true = abs([Paths_true.alpha]);
alpha_true_dB = 20*log10(alpha_true / max(alpha_true));

alpha_est = abs(ParaEst.Alpha_MPCE(1,:));
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


end
