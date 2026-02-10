
function alpha_hat = Z_Corr_mx_MPCE(ParaEst, Le, Dataxl, ChPara)


Pa.AoD_l = ParaEst.thetaTX_MPCE(Le);
Pa.AoA_l = ParaEst.thetaRX_MPCE(Le);
Pa.tau_l = ParaEst.DelayE(Le);
        
% Constructing Angle Time Space of MPC
A_Beamsteering_delay = Construction_ofMPC(Pa, ChPara);

alpha_hat = sum(conj(A_Beamsteering_delay(:)) .* Dataxl(:)) / (sum(abs(A_Beamsteering_delay(:)).^2));  % complex gain estimate     


end

