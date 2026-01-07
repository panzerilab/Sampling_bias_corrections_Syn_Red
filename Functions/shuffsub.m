function lattice_corr = shuffsub(input_data, nshuff, pid)
% SHUFFSUB  Shuffle–substitution control for PID lattice.
% input_data: V x N matrix of observations
%   - Discrete: integer labels (1..K) per variable
%   - Gaussian: real-valued samples
% nshuff:     number of shuffles
% pid:        pid_lattice object (uses pid.is_gaussian, pid.nsources)
%
% Returns:
%   lattice_corr = lattice - lattice_shuff   (size: nAtoms x nshuff)

    if nargin < 3
        error('Provide the pid_lattice object as the 3rd input.');
    end

    % N x V (rows = trials, cols = variables [sources..., target])
    obs = input_data.'; 
    [nTrials, total_vars] = size(obs);

    % Identify target column(s)
    if isprop(pid, 'target_dims') && ~isempty(pid.target_dims)
        tgt = pid.target_dims(:).';           % allow vector target
    else
        tgt = pid.nsources + 1;               % default: last column
    end
    if any(tgt < 1 | tgt > total_vars)
        error('Target column(s) out of range for provided data.');
    end

    if pid.is_gaussian
        % ---------- Gaussian path ----------
        % Lattice on observed covariance
        Sigma = cov(obs, 1);
        Sigma = (Sigma + Sigma.')/2;
        lattice = pid.calculate_latvals(Sigma);

        L = numel(lattice);
        lattice_shuff = zeros(L, nshuff);

        for ns = 1:nshuff
            shuf = obs;
            % Shuffle target rows jointly (preserve target's own structure)
            rp = randperm(nTrials);
            shuf(:, tgt) = obs(rp, tgt);

            Sigma_res = cov(shuf, 1);
            Sigma_res = (Sigma_res + Sigma_res.')/2;
            lattice_shuff(:, ns) = pid.calculate_latvals(Sigma_res);
        end

    else
        % ---------- Discrete path ----------
        % Ensure 1-based labels
        mins = min(obs, [], 1);
        if any(mins < 1)
            obs = bsxfun(@minus, obs, mins) + 1;
        end

        % Size of joint table per variable
        nb_vec = max(obs, [], 1);

        % Observed PMF
        counts = accumarray(obs, 1, nb_vec);
        p = counts / sum(counts(:));
        lattice = pid.calculate_latvals(p);

        L = numel(lattice);
        lattice_shuff = zeros(L, nshuff);

        for ns = 1:nshuff
            shuf = obs;
            % Shuffle only the target column (rows), not columns!
            rp = randperm(nTrials);
            shuf(:, tgt) = obs(rp, tgt);

            % Rebuild PMF from shuffled table
            counts_res = accumarray(shuf, 1, nb_vec);
            p_res = counts_res / sum(counts_res(:));
            lattice_shuff(:, ns) = pid.calculate_latvals(p_res);
        end
    end

    % Per-shuffle difference (same convention you used)
    lattice_corr = lattice - mean(lattice_shuff,2)';
end
