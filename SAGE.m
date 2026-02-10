function Parameter_est = SAGE(Measured_MPC, No_iteration, ChPara)

DataTimeDomain = Measured_MPC;

% ---- Initialization (SIC) ----
thetaTX_MPCE = zeros(1, ChPara.N);
thetaRX_MPCE = zeros(1, ChPara.N);
DelayE       = zeros(1, ChPara.N);
mx_MPCE      = zeros(1, ChPara.N);

DataTimeDomain_MPC = Measured_MPC;       % Ptx x Prx x Ndelay,  residual for SIC

for Li = 1:ChPara.N
    Param = make_initial_guess(DataTimeDomain_MPC, ChPara);

    thetaTX_MPCE(Li) = Param.AoD_l;
    thetaRX_MPCE(Li) = Param.AoA_l;
    DelayE(Li)       = Param.tau_l;
    mx_MPCE(Li)      = Param.alpha_l;
    


    % Subtract estimated MPC (SIC)
%     Sl_0  = Param.alphaMPC0*Construction_ofMPC(Param, ChPara);
%      DataTimeDomain_MPC = DataTimeDomain_MPC - Sl_0;
    DataTimeDomain_MPC = DataTimeDomain_MPC -   Param.alpha_l*Construction_ofMPC(Param, ChPara);
    


end

ParaEst = struct();
ParaEst.thetaTX_MPCE = thetaTX_MPCE;
ParaEst.thetaRX_MPCE = thetaRX_MPCE;
ParaEst.DelayE       = DelayE;
ParaEst.mx_MPCE      = mx_MPCE;

disp('--- Initialization (MPC Parameters) ---');

alpha_hat = abs(ParaEst.mx_MPCE);
% Due to steering matrix normalization, 
% the estimated path gains are reported in relative form, 
% normalized with respect to the strongest MPC
alpha_rel = alpha_hat / max(alpha_hat);



for i = 1:ChPara.N
    fprintf('MPC %d: AoD = %6.1f deg, AoA = %6.1f deg, Delay = %6.2f ns, |α|_rel = %.3f\n', ...
        i, ...
        ParaEst.thetaTX_MPCE(i)*180/pi, ...
        ParaEst.thetaRX_MPCE(i)*180/pi, ...
        ParaEst.DelayE(i)*ChPara.Ts*1e9, ...   % convert to ns
        alpha_rel(i));
end


% ---- EM / SAGE iterations ----

accept_update = @(Znew, Zold) (Znew > Zold);

% simple convergence thresholds
th_Ang_deg = 1;   % degrees
th_Delay   = 0.01;  % samples

for n = 1:No_iteration

    % save previous parameters
    thetaTX_prev = ParaEst.thetaTX_MPCE(:);
    thetaRX_prev = ParaEst.thetaRX_MPCE(:);
    delay_prev   = ParaEst.DelayE(:);

    for Le = 1:ChPara.N

        % E-step
        Dataxl = ExpectationXL(DataTimeDomain, ParaEst, Le, ChPara);

        % AoD update
        [thetaTX_new, ZTX_new, ZTX_old] = Z_Corr_thetaTX(ParaEst, Le, Dataxl, ChPara);
        if accept_update(ZTX_new, ZTX_old)
            ParaEst.thetaTX_MPCE(Le) = thetaTX_new;
        end

        % AoA update
        [thetaRX_new, ZRX_new, ZRX_old] = Z_Corr_thetaRX(ParaEst, Le, Dataxl, ChPara);
        if accept_update(ZRX_new, ZRX_old)
            ParaEst.thetaRX_MPCE(Le) = thetaRX_new;
        end

        % Delay update
        [d_new, ZD_new, ZD_old] = Z_Corr_Delay(ParaEst, Le, Dataxl, ChPara);
        if accept_update(ZD_new, ZD_old)
            ParaEst.DelayE(Le) = d_new;
        end

        % amplitude
        ParaEst.mx_MPCE(Le) = Z_Corr_mx_MPCE(ParaEst, Le, Dataxl, ChPara);

    end

    % ---- display ----
    disp(['--- Iteration ', num2str(n), ' (MPC Parameters) ---']);

    alpha_hat = abs(ParaEst.mx_MPCE);
    alpha_rel = alpha_hat / max(alpha_hat);

    for i = 1:ChPara.N
        fprintf('MPC %d: AoD = %6.1f deg, AoA = %6.1f deg, Delay = %6.2f ns, |α|_rel = %.3f\n', ...
            i, ParaEst.thetaTX_MPCE(i)*180/pi, ...
            ParaEst.thetaRX_MPCE(i)*180/pi, ...
            ParaEst.DelayE(i)*ChPara.Ts*1e9, ...
            alpha_rel(i));
    end

    % ---- simple convergence check ----
    max_dAoD = max(abs(ParaEst.thetaTX_MPCE(:) - thetaTX_prev(:))) * 180/pi;
    max_dAoA = max(abs(ParaEst.thetaRX_MPCE(:) - thetaRX_prev(:))) * 180/pi;
    max_dDel = max(abs(ParaEst.DelayE(:)       - delay_prev(:)));


    fprintf('Max ΔAoD = %.3f deg, Max ΔAoA = %.3f deg, Max ΔDelay = %.3f samples\n', ...
        max_dAoD, max_dAoA, max_dDel);

    if (max_dAoD < th_Ang_deg) && (max_dAoA < th_Ang_deg) && (max_dDel < th_Delay)
        disp(['Converged at iteration ', num2str(n)]);
        break;
    end

end


Parameter_est = ParaEst;

end
