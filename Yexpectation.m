function DataTimeDomain_ = Yexpectation(DataTimeDomain_, Pa, ChPara)

thetaRX = ChPara.Ang_steer;
thetaTX = ChPara.Ang_steer;

S = zeros(length(thetaTX), length(thetaRX), length(ChPara.f));

TX_Chain = struct();
RX_Chain = struct();

% TX steering
for z = 1:length(ChPara.Ang_steer)
    frx = exp(-1i*(ChPara.n-1)*2*pi*ChPara.d*(cos(thetaTX(z)) - cos(Pa.thetaTX_MPC))/ChPara.lambda);
    TX_Chain.RF1(z) = sum(frx);
end
RF_vector_Tx = [TX_Chain.RF1];

% RX steering
for z = 1:length(ChPara.Ang_steer)
    frx = exp(-1i*(ChPara.n-1)*2*pi*ChPara.d*(cos(thetaRX(z)) - cos(Pa.thetaRX_MPC))/ChPara.lambda);
    RX_Chain.RF1(z) = sum(frx);
end
RF_vector_Rx = [RX_Chain.RF1];

% A(thetaT,thetaR)
AmRFnRF = conj(RF_vector_Tx.') * RF_vector_Rx;

% normalize A (recommended for measurement beamsteering)
% maxA = max(abs(AmRFnRF(:)));
% AmRFnRF = AmRFnRF / maxA;


% Delay impulse in delay domain (Pa.Delay_MPC is in samples)
Delay_domain = ifft(exp(-1i * 2*pi * (ChPara.f) * Pa.Delay_MPC * ChPara.Ts));

% Build estimated MPC in delay domain: alpha * A * delaytap
for k = 1:length(Delay_domain)
    S(:,:,k) = Pa.mxMPC * (AmRFnRF * Delay_domain(k));
end

% Subtract (path cancellation)
DataTimeDomain_ = DataTimeDomain_ - S;

end
