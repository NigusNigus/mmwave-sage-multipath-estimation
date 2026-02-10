function Xl = ExpectationXL(DataTimeDomain_, ParaEst, Le, ChPara)

% SAGE E-step: construction of path-specific residual X_l
%   X_l = X - sum_{l' ~= l} alpha_l' * A(thetaT_l', thetaR_l') * delay_domain(tau_l')
% This function removes the contribution of all paths except the l-th one,
% producing the residual Xl used in the M-step to update path l.

thetaRX = ChPara.Ang_steer;
thetaTX = ChPara.Ang_steer;

A_Beamsteering_delay = zeros(length(thetaTX), length(thetaRX), length(ChPara.f));
A_Beamsteering_delay_MPCs = zeros(size(A_Beamsteering_delay));

for l = 1:ChPara.N

    if l ~= Le

        Pa.AoD_l = ParaEst.thetaTX_MPCE(l);
        Pa.AoA_l = ParaEst.thetaRX_MPCE(l);
        Pa.tau_l = ParaEst.DelayE(l);
        Pa.alpha_l =ParaEst.mx_MPCE(l);
        
        % Constructing Angle Time Space of MPC
        A_Beamsteering_delay = Construction_ofMPC(Pa, ChPara);

        
        A_Beamsteering_delay_MPCs = A_Beamsteering_delay_MPCs +  Pa.alpha_l*A_Beamsteering_delay;



    end
    
end


% Residual for path Le
Xl = DataTimeDomain_ - A_Beamsteering_delay_MPCs;


end
