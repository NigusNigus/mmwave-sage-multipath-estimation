function Param = make_initial_guess(DataTimeDomain, ChPara)

Param = struct();
DataTimeDomain_ = DataTimeDomain;

% Coarse AoD/AoA from global peak (from Tx/Rx Beamsteerings Angles)
[~, max_index] = max(abs(DataTimeDomain_(:)));
[aod_index, aoa_index, ~] = ind2sub(size(DataTimeDomain_), max_index);

Param.AoD_l = ChPara.Ang_steer(aod_index);
Param.AoA_l = ChPara.Ang_steer(aoa_index);

%  Delay refinement at that (AoD,AoA)
DelayGrid = 0:0.01:10; % delays to evaluate, samples,  used 10 to save the iteration, should be 100
f = ChPara.f(:);


Yt= DataTimeDomain_(aod_index, aoa_index, :);
Yt = Yt(:);

Dela_corr_Resp_t = zeros(1, numel(DelayGrid));
for k = 1:numel(DelayGrid)
    Delay = DelayGrid(k); % samples
    dt = ifft(exp(-1i*2*pi*f*Delay*ChPara.Ts)); 
    Dela_corr_Resp_t(k) = dt'*Yt;
end


[~, maxIdx] = max(abs(Dela_corr_Resp_t));
Delay_hat = DelayGrid(maxIdx);
Param.tau_l = Delay_hat;   % samples (tau = Delay_hat*Ts)


%Intilization of alpha_MPC
A_beamsteering_delay_est = Construction_ofMPC(Param, ChPara);
Param.alpha_l = sum(conj(A_beamsteering_delay_est(:)) .* DataTimeDomain_(:)) / (sum(abs(A_beamsteering_delay_est(:)).^2));


end
