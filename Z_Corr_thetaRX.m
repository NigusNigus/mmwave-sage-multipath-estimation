function [Para, Znew, Zold] = Z_Corr_thetaRX(ParaEst, Le, Dataxl,  ChPara)

thetaRX_MPCE_K = -45*pi/180:1*pi/180:45*pi/180;


Zvector = zeros(1, numel(thetaRX_MPCE_K));


Pa.AoD_l = ParaEst.thetaTX_MPCE(Le);
Pa.tau_l = ParaEst.DelayE(Le);

% Zold at current AoA
Pa.AoA_l = ParaEst.thetaRX_MPCE(Le);
A_old = Construction_ofMPC(Pa, ChPara); %Beamsteering-Response for old MPC parameters
Zold = abs(sum(conj(A_old(:)) .* Dataxl(:)));

% Search for best AoA
for i_angle = 1:numel(thetaRX_MPCE_K)
    
    Pa.AoA_l = thetaRX_MPCE_K(i_angle);

    % Constructing Angle Time Space of MPC
    A_Beamsteering_delay = Construction_ofMPC(Pa, ChPara);

    Zvector(i_angle) = abs(sum(conj(A_Beamsteering_delay(:)) .* Dataxl(:)));

end

[Znew, idx] = max(Zvector);
Para = thetaRX_MPCE_K(idx);

end


