function Paths = build_paths(Tx, Rx_list, Ref_list, Segments, Meta)
    NRx  = size(Rx_list,1);
    NRef = size(Ref_list,1);
    

    angfun = @(p1,p2) atan2(p2(1)-p1(1), p2(2)-p1(2));
    c = 3e8;
    fc = 28e9;   % 28 GHz

    Paths = struct('rxIdx',{},'type',{},'refIdx',{}, ...
                   'segTxRef',{},'segRefRx',{},'segLOS',{}, ...
                   'segLast',{},'thetaTX',{},'thetaRX',{}, ...
                   'tau',{},'alpha',{},'phi',{});

    idx = 0;

    % LOS paths
    for j = 1:NRx
        idx = idx+1;
        segLOS = find(strcmp({Meta.type},'TxRx') & [Meta.rxIdx]==j);
        if numel(segLOS) ~= 1
            error('LOS segment for Rx%d not found/unique', j);
        end
        Paths(idx).rxIdx    = j;
        Paths(idx).type     = 'LOS';
        Paths(idx).refIdx   = NaN;
        Paths(idx).segTxRef = NaN;
        Paths(idx).segRefRx = NaN;
        Paths(idx).segLOS   = segLOS;
        Paths(idx).segLast  = segLOS;

        th = angfun(Tx, Rx_list(j,:));
        Paths(idx).thetaTX = th;
        Paths(idx).thetaRX = th;

        L = norm(Rx_list(j,:)-Tx);
        Paths(idx).tau   = L/c;
        phi = -2*pi*fc*Paths(idx).tau;
        Paths(idx).phi = mod(phi, 2*pi);
        Paths(idx).alpha = (1/L) * exp(1j*phi);
      
    end

    % Reflected paths Tx->Ref->Rx
    for j = 1:NRx
        for r = 1:NRef
            idx = idx+1;
            segTxRef = find(strcmp({Meta.type},'TxRef') & [Meta.refIdx]==r);
            segRefRx = find(strcmp({Meta.type},'RefRx') & ...
                            [Meta.refIdx]==r & [Meta.rxIdx]==j);
            if numel(segTxRef)~=1 || numel(segRefRx)~=1
                error('Segments for Tx->Ref%d->Rx%d not found/unique', r, j);
            end
            Paths(idx).rxIdx    = j;
            Paths(idx).type     = 'REF';
            Paths(idx).refIdx   = r;
            Paths(idx).segTxRef = segTxRef;
            Paths(idx).segRefRx = segRefRx;
            Paths(idx).segLOS   = NaN;
            Paths(idx).segLast  = segRefRx;

            thTX = angfun(Tx, Ref_list(r,:));
            thRX = angfun(Ref_list(r,:), Rx_list(j,:));
            Paths(idx).thetaTX = thTX;
            Paths(idx).thetaRX = thRX;

            L = norm(Ref_list(r,:)-Tx) + norm(Rx_list(j,:)-Ref_list(r,:));
            Paths(idx).tau   = L/c;
            Paths(idx).tau   = L/c;
            phi = -2*pi*fc*Paths(idx).tau;
            Paths(idx).phi = mod(phi, 2*pi);
            Paths(idx).alpha = (1/L) * exp(1j*phi);
            
        end
    end
end