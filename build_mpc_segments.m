function [Segments, Names, Meta] = build_mpc_segments(Tx, Rx_list, Ref_list)
    NRx  = size(Rx_list,1);
    NRef = size(Ref_list,1);

    Segments = [];
    Names    = {};
    Meta     = struct('type',{},'refIdx',{},'rxIdx',{});

    % Tx->Ref
    for r = 1:NRef
        Segments(end+1,:) = [Tx Ref_list(r,:)]; %#ok<AGROW>
        Names{end+1} = sprintf('Tx-Ref%d',r);
        Meta(end+1).type   = 'TxRef';
        Meta(end).refIdx   = r;
        Meta(end).rxIdx    = NaN;
    end

    % Ref->Rx
    for r = 1:NRef
        for j = 1:NRx
            Segments(end+1,:) = [Ref_list(r,:) Rx_list(j,:)]; %#ok<AGROW>
            Names{end+1} = sprintf('Ref%d-Rx%d',r,j);
            Meta(end+1).type   = 'RefRx';
            Meta(end).refIdx   = r;
            Meta(end).rxIdx    = j;
        end
    end

    % Tx->Rx (LOS)
    for j = 1:NRx
        Segments(end+1,:) = [Tx Rx_list(j,:)]; %#ok<AGROW>
        Names{end+1} = sprintf('Tx-Rx%d',j);
        Meta(end+1).type   = 'TxRx';
        Meta(end).refIdx   = NaN;
        Meta(end).rxIdx    = j;
    end
end