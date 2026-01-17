function [pid_v] = GroundTruth(edges_spikes1uc, edges_spikes2uc, rate1, rate2, rateSh, nRbins, nStimuli, numX, pidMeasure)

p_s = 1/nStimuli;
p_x1x2_s = zeros(nStimuli, numX, numX);
p_x1x2s = zeros(nStimuli, numX, numX);

totalComb = nStimuli * numX * numX;
results = zeros(totalComb, 1);
indexMap = zeros(totalComb, 3); % Store [s, x1+1, x2+1]

% Precompute combinations
[S, X1, X2] = ndgrid(1:nStimuli, 0:numX-1, 0:numX-1);
S = S(:); X1 = X1(:); X2 = X2(:);

tmp_vals = zeros(totalComb, 1);
switch pidMeasure
    case 'I_BROJA'
        PIDfunction = @pidBROJA;
    case 'I_MMI'
        PIDfunction = @pidimmi;
    case 'I_min'
        PIDfunction = @pidimin;
end
parfor idx = 1:totalComb
    s = S(idx);
    x1 = X1(idx);
    x2 = X2(idx);

    p_val = 0;
    for xsh = 0:min(x1, x2)
        p_rSh = poisspdf(xsh, rateSh(s));
        X1_tmp = x1 - xsh;
        X2_tmp = x2 - xsh;
        p1_tmp = poisspdf(X1_tmp, rate1(s));
        p2_tmp = poisspdf(X2_tmp, rate2(s));
        p_val = p_val + (p_rSh * p1_tmp * p2_tmp);
    end

    tmp_vals(idx) = p_val;
    indexMap(idx, :) = [s, x1 + 1, x2 + 1];
end

% Accumulate results back into arrays
for idx = 1:totalComb
    s = indexMap(idx, 1);
    x1_idx = indexMap(idx, 2);
    x2_idx = indexMap(idx, 3);
    p_x1x2_s(s, x1_idx, x2_idx) = tmp_vals(idx);
    p_x1x2s(s, x1_idx, x2_idx) = tmp_vals(idx) * p_s;
end

% ---- Binning and PID (unchanged) ----
if length(nRbins) > 1
    for binIDX = 1:length(nRbins)
        rbinNum = nRbins(binIDX);
        edges_spikes1uc_tmp = edges_spikes1uc{binIDX};
        edges_spikes2uc_tmp = edges_spikes2uc{binIDX};

        if edges_spikes1uc_tmp(end-1) < numX && edges_spikes2uc_tmp(end-1) < numX
            edges_spikes1uc_tmp(end) = numX;
            edges_spikes2uc_tmp(end) = numX;
        end

        p_x1x2_s_bins = zeros(nStimuli, rbinNum, rbinNum);
        p_x1x2s_bins = zeros(nStimuli, rbinNum, rbinNum);
        p_x1x2 = zeros(nStimuli, rbinNum, rbinNum);

        for s = 1:nStimuli
            for rB1 = 1:rbinNum
                for rB2 = 1:rbinNum
                    range1 = (edges_spikes1uc_tmp(rB1)+1):(edges_spikes1uc_tmp(rB1+1));
                    range2 = (edges_spikes2uc_tmp(rB2)+1):(edges_spikes2uc_tmp(rB2+1));
                    p_x1x2_s_bins(s, rB1, rB2) = sum(sum(p_x1x2_s(s, range1, range2)));
                    p_x1x2s_bins(s, rB1, rB2) = sum(sum(p_x1x2s(s, range1, range2)));
                    p_x1x2(s, rB1, rB2) = sum(sum(p_x1x2s(s, range1, range2)));
                end
            end
        end

        P = permute(p_x1x2s_bins, [2 3 1]);
        % Compute P(S)
        P_S = squeeze(sum(P, [1, 2])); % Summing over X1 and X2
        X1_size =size(P,1);
        X2_size =size(P,2);
        S_size  =size(P,3);
        % Compute P(X1 | S) and P(X2 | S)
        P_X1_given_S = zeros(X1_size, S_size);
        P_X2_given_S = zeros(X2_size, S_size);

        for s = 1:S_size
            P_X1_given_S(:, s) = sum(P(:, :, s), 2) ./ P_S(s); % Conditional P(X1 | S)
            P_X2_given_S(:, s) = sum(P(:, :, s), 1) ./ P_S(s); % Conditional P(X2 | S)
        end

        % Compute Pind(X1, X2 | S)
        Pind_X1_X2_given_S = zeros(X1_size, X2_size, S_size);
        for s = 1:S_size
            for x1 = 1:X1_size
                for x2 = 1:X2_size
                    Pind_X1_X2_given_S(x1, x2, s) = P_X1_given_S(x1, s) * P_X2_given_S(x2, s);
                end
            end
        end

        % Convert to joint probability using P(S)
        Pind = zeros(X1_size, X2_size, S_size);
        for s = 1:S_size
            Pind(:, :, s) = Pind_X1_X2_given_S(:, :, s) * P_S(s);
        end

        % Normalize Pind to ensure it is a valid probability distribution
        Pind = Pind ./ sum(Pind, 'all');
        Pind =  permute(Pind, [3 1 2]);
        Pind_v_tmp = mutualInformationXYZ(Pind);


        % pid_v_tmp = pidBROJA(p_x1x2);
        % pid_v{binIDX} = pid_v_tmp;
        pid_v_tmp = PIDfunction(p_x1x2);
        % qdist_v_tmp = mutualInformationXYZ(qdist_tmp);
        pid_v_tmp = [pid_v_tmp, Pind_v_tmp];
        pid_v{binIDX} = pid_v_tmp;
    end
else
    if edges_spikes1uc(end-1) < numX && edges_spikes2uc(end-1) < numX
        edges_spikes1uc(end) = numX;
        edges_spikes2uc(end) = numX;
    end

    p_sx1x2_bins = zeros(nStimuli, nRbins, nRbins);
    p_x1x2 = zeros(nStimuli, nRbins, nRbins);

    for s = 1:nStimuli
        for rB1 = 1:nRbins
            for rB2 = 1:nRbins
                range1 = (edges_spikes1uc(rB1)+1):(edges_spikes1uc(rB1+1));
                range2 = (edges_spikes2uc(rB2)+1):(edges_spikes2uc(rB2+1));
                p_sx1x2_bins(s, rB1, rB2) = sum(sum(p_x1x2s(s, range1, range2)));
                p_x1x2(s, rB1, rB2)       = sum(sum(p_x1x2s(s, range1, range2)));
            end
        end
    end

    switch pidMeasure
        case 'I_BROJA'
            PIDfunction = @pidBROJA;
        case 'I_MMI'
            PIDfunction = @pidimmi;
        case 'I_min'
            PIDfunction = @pidimin;
    end

    P = permute(p_sx1x2_bins, [2 3 1]);
    % Compute P(S)
    P_S = squeeze(sum(P, [1, 2])); % Summing over X1 and X2
    X1_size =size(P,1);
    X2_size =size(P,2);
    S_size  =size(P,3);
    % Compute P(X1 | S) and P(X2 | S)
    P_X1_given_S = zeros(X1_size, S_size);
    P_X2_given_S = zeros(X2_size, S_size);

    for s = 1:S_size
        P_X1_given_S(:, s) = sum(P(:, :, s), 2) ./ P_S(s); % Conditional P(X1 | S)
        P_X2_given_S(:, s) = sum(P(:, :, s), 1) ./ P_S(s); % Conditional P(X2 | S)
    end

    % Compute Pind(X1, X2 | S)
    Pind_X1_X2_given_S = zeros(X1_size, X2_size, S_size);
    for s = 1:S_size
        for x1 = 1:X1_size
            for x2 = 1:X2_size
                Pind_X1_X2_given_S(x1, x2, s) = P_X1_given_S(x1, s) * P_X2_given_S(x2, s);
            end
        end
    end

    % Convert to joint probability using P(S)
    Pind = zeros(X1_size, X2_size, S_size);
    for s = 1:S_size
        Pind(:, :, s) = Pind_X1_X2_given_S(:, :, s) * P_S(s);
    end

    % Normalize Pind to ensure it is a valid probability distribution
    Pind = Pind ./ sum(Pind, 'all');
    Pind =  permute(Pind, [3 1 2]);
    Pind_v = mutualInformationXYZ(Pind);

    pid_v = PIDfunction(p_x1x2);
    % qdist_v = mutualInformationXYZ(qdist);
    pid_v = [pid_v, Pind_v];
end
end
