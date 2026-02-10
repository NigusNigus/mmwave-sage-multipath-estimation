function [Para, Znew, Zold] = Z_Corr_Delay(ParaEst, Le, Dataxl, ChPara)

DelayE_MPCE_K = 0:0.01:10;  % The maximum is 100

Zvector = zeros(1, numel(DelayE_MPCE_K));

Pa.AoD_l = ParaEst.thetaTX_MPCE(Le);
Pa.AoA_l = ParaEst.thetaRX_MPCE(Le);

% --- Zold at current delay
Pa.tau_l = ParaEst.DelayE(Le);
A_old = Construction_ofMPC(Pa, ChPara);  %Beamsteering-Response for old MPC parameters
Zold = abs(sum(conj(A_old(:)) .* Dataxl(:)));
 

% Search for best delay
for i_delay = 1:numel(DelayE_MPCE_K)

    Pa.tau_l =DelayE_MPCE_K(i_delay); 

    % Constructing Angle Time Space of MPC
    A_Beamsteering_delay = Construction_ofMPC(Pa, ChPara);

    Zvector(i_delay) = abs(sum(conj(A_Beamsteering_delay(:)) .* Dataxl(:)));

end

[Znew, idx] = max(Zvector);

% NOTE: Para is the delay estimate in samples.
% Physical propagation delay (seconds) is Para * ChPara.Ts.
Para = DelayE_MPCE_K(idx);

end

