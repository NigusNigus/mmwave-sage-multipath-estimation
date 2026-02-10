function A_beamsteering_delay = Construction_ofMPC(Pa, ChPara)

thetaRX = ChPara.Ang_steer;
thetaTX = ChPara.Ang_steer;

A_beamsteering_delay = zeros(length(thetaTX), length(thetaRX), length(ChPara.f));

TX_Chain = struct();
RX_Chain = struct();

% TX steering
for z = 1:length(ChPara.Ang_steer)
    frx = exp(-1i*(ChPara.n)*2*pi*ChPara.d*(sin(thetaTX(z)) - sin(Pa.AoD_l))/ChPara.lambda);
    TX_Chain.RF1(z) = sum(frx);
end
RF_vector_Tx = [TX_Chain.RF1];

% RX steering
for z = 1:length(ChPara.Ang_steer)
    frx = exp(-1i*(ChPara.n)*2*pi*ChPara.d*(sin(thetaRX(z)) - sin(Pa.AoA_l))/ChPara.lambda);
    RX_Chain.RF1(z) = sum(frx);
end
RF_vector_Rx = [RX_Chain.RF1];

% A(thetaT,thetaR)
A_beamsteering = conj(RF_vector_Tx.') * RF_vector_Rx;


% Delay impulse in delay domain (Pa.Delay_MPC is in samples)
Delay_domain = ifft(exp(-1i * 2*pi * (ChPara.f) * Pa.tau_l * ChPara.Ts));


% Build estimated MPC in delay domain: space–angle–delay steering matrix 
for k = 1:length(Delay_domain)
    A_beamsteering_delay(:,:,k) =(A_beamsteering * Delay_domain(k));
end


end