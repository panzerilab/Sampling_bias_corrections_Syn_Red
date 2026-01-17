function lattice_corr = resample(input_data, nshuff, pid)
% RESAMPLE  Bias-correct PID lattice via parametric bootstrap.
% input_data: V x N matrix of observations
%   - Discrete case: integer state labels (1..K) for each variable
%   - Gaussian case: real-valued samples (sources + target)
% nshuff:     number of bootstrap resamples
% pid:        pid_lattice object (reads pid.is_gaussian, nsources, etc.)
%
% Returns:
%   lattice_corr  = 2*lattice - mean(lattice_shuff,2)

    if nargin < 3
        error('Provide the pid_lattice object as the 3rd input.');
    end

    obs = input_data.';                     % N x V
    [nTrials, total_vars] = size(obs);

    % --- Evaluate on observed data ---
    if pid.is_gaussian
        % Mean/cov of continuous data (cov(...,1) -> 1/N normalization)
        mu = mean(obs, 1);
        Sigma = cov(obs, 1);
        Sigma = (Sigma + Sigma.')/2;       % symmetrize

        lattice = pid.calculate_latvals(Sigma);
        L = numel(lattice);
        lattice_shuff = zeros(L, nshuff);

        % Factor once (with jitter) for stable MVN sampling
        R = chol_with_jitter(Sigma);

        for ns = 1:nshuff
            Z = randn(nTrials, total_vars) * R' + mu;  % MVN(nTrials, mu, Sigma)
            Sigma_res = cov(Z, 1);
            Sigma_res = (Sigma_res + Sigma_res.')/2;
            lattice_shuff(:, ns) = pid.calculate_latvals(Sigma_res);
        end

    else
        % ------ Discrete (multinomial) path ------
        % Shift to 1-based if needed
        mins = min(obs, [], 1);
        if any(mins < 1)
            obs = bsxfun(@minus, obs, mins) + 1;
        end

        % Bins per variable
        nb_vec = max(obs, [], 1);

        % Joint counts and PMF
        counts = accumarray(obs, 1, nb_vec);
        p = counts / sum(counts(:));

        lattice = pid.calculate_latvals(p);
        L = numel(lattice);
        lattice_shuff = zeros(L, nshuff);

        % Flatten PMF for sampling
        pmf = p(:);
        pmf = pmf / sum(pmf);
        num_states = numel(pmf);

        for ns = 1:nshuff
            % Sample joint states from fitted PMF
            idx = randsample(num_states, nTrials, true, pmf);

            % Convert to per-variable subscripts (N x V)
            subs = cell(1, total_vars);
            [subs{:}] = ind2sub(nb_vec, idx);
            resampled_obs = cell2mat(subs);

            % Re-estimate PMF for the resample
            counts_res = accumarray(resampled_obs, 1, nb_vec);
            p_res = counts_res / sum(counts_res(:));

            lattice_shuff(:, ns) = pid.calculate_latvals(p_res);
        end
    end

    % Bias-corrected estimate
    lattice_corr = 2 * lattice - mean(lattice_shuff, 2)';
end

% ---------- helpers ----------
function R = chol_with_jitter(S)
    % Cholesky with escalating jitter to ensure SPD
    [R, p] = chol(S, 'lower');
    if p == 0, return; end
    jitter = 1e-10;
    I = eye(size(S));
    for k = 1:8
        [R, p] = chol(S + jitter*I, 'lower');
        if p == 0, return; end
        jitter = jitter * 10;
    end
    error('Covariance not SPD even after jittering.');
end
