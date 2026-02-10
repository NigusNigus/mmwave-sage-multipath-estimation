function [Para, Znew, Zold] = Z_Corr_thetaTX(ParaEst, Le, Dataxl, ChPara)

thetaTX_MPCE_K = -50*pi/180:1*pi/180:50*pi/180;

Zvector = zeros(1, numel(thetaTX_MPCE_K));

Pa.AoA_l = ParaEst.thetaRX_MPCE(Le);
Pa.tau_l = ParaEst.DelayE(Le);

% Zold at current AoD
Pa.AoD_l = ParaEst.thetaTX_MPCE(Le);
A_old = Construction_ofMPC(Pa, ChPara); %Beamsteering-Response for old MPC parameters
Zold = abs(sum(conj(A_old(:)) .* Dataxl(:)))/ (sum(abs(A_old(:)).^2));

% Search for best AoD
for i_angle = 1:numel(thetaTX_MPCE_K)

    Pa.AoD_l = thetaTX_MPCE_K(i_angle);

    % Constructing Angle Time Space of MPC
    A_Beamsteering_delay = Construction_ofMPC(Pa, ChPara);

    Zvector(i_angle) = abs(sum(conj(A_Beamsteering_delay(:)) .* Dataxl(:))) / (sum(abs(A_Beamsteering_delay(:)).^2));

end

[Znew, idx] = max(Zvector);
Para = thetaTX_MPCE_K(idx);


end
