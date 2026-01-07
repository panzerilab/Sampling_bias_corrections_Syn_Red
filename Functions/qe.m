function lattice_corr = qe(input_data, nshuff, pid)
% QE  Quadratic-extrapolation bias correction for PID lattice.
% input_data: V x N observations
%   - Discrete: integer labels (1..K) per variable
%   - Gaussian: real-valued samples (sources then target)
% nshuff:     number of random partitions
% pid:        pid_lattice object (reads pid.is_gaussian, pid.nsources, etc.)
%
% Returns:
%   lattice_corr  : bias-corrected lattice atoms via quadratic extrapolation
%                   over 1/N using full, half, and quarter sample sizes.

    if nargin < 3
        error('Provide the pid_lattice object as the 3rd input.');
    end

    % Work in N x V (rows = trials)
    obs = input_data.'; 
    [nTrials, total_vars] = size(obs);

    if nTrials < 4
        error('QE requires at least 4 trials (got %d).', nTrials);
    end

    % ---------- Evaluate lattice on full sample ----------
    if pid.is_gaussian
        Sigma = cov(obs, 1);
        Sigma = (Sigma + Sigma.')/2;
        lattice = pid.calculate_latvals(Sigma);

        L = numel(lattice);
        lattice2_shuff = zeros(L, 2, nshuff);  % halves
        lattice4_shuff = zeros(L, 4, nshuff);  % quarters

        for ns = 1:nshuff
            perm = randperm(nTrials);
            n2 = floor(nTrials/2);
            n4 = floor(nTrials/4);

            h1 = perm(1:n2);
            h2 = perm(n2+1:2*n2);

            q1 = perm(1:n4);
            q2 = perm(n4+1:2*n4);
            q3 = perm(2*n4+1:3*n4);
            q4 = perm(3*n4+1:4*n4);

            lattice2_shuff(:,1,ns) = lat_from_rows_gauss(obs, h1, pid);
            lattice2_shuff(:,2,ns) = lat_from_rows_gauss(obs, h2, pid);

            lattice4_shuff(:,1,ns) = lat_from_rows_gauss(obs, q1, pid);
            lattice4_shuff(:,2,ns) = lat_from_rows_gauss(obs, q2, pid);
            lattice4_shuff(:,3,ns) = lat_from_rows_gauss(obs, q3, pid);
            lattice4_shuff(:,4,ns) = lat_from_rows_gauss(obs, q4, pid);
        end

    else
        % -------- Discrete path: build PMF via accumarray --------
        % Ensure 1-based labels
        mins = min(obs, [], 1);
        if any(mins < 1)
            obs = bsxfun(@minus, obs, mins) + 1;
        end
        nb_vec = max(obs, [], 1);

        counts = accumarray(obs, 1, nb_vec);
        p = counts / sum(counts(:));
        lattice = pid.calculate_latvals(p);

        L = numel(lattice);
        lattice2_shuff = zeros(L, 2, nshuff);
        lattice4_shuff = zeros(L, 4, nshuff);

        for ns = 1:nshuff
            perm = randperm(nTrials);
            n2 = floor(nTrials/2);
            n4 = floor(nTrials/4);

            h1 = perm(1:n2);
            h2 = perm(n2+1:2*n2);

            q1 = perm(1:n4);
            q2 = perm(n4+1:2*n4);
            q3 = perm(2*n4+1:3*n4);
            q4 = perm(3*n4+1:4*n4);

            lattice2_shuff(:,1,ns) = lat_from_rows_disc(obs, h1, nb_vec, pid);
            lattice2_shuff(:,2,ns) = lat_from_rows_disc(obs, h2, nb_vec, pid);

            lattice4_shuff(:,1,ns) = lat_from_rows_disc(obs, q1, nb_vec, pid);
            lattice4_shuff(:,2,ns) = lat_from_rows_disc(obs, q2, nb_vec, pid);
            lattice4_shuff(:,3,ns) = lat_from_rows_disc(obs, q3, nb_vec, pid);
            lattice4_shuff(:,4,ns) = lat_from_rows_disc(obs, q4, nb_vec, pid);
        end
    end

    % ---------- Quadratic extrapolation over x = [1/N, 2/N, 4/N] ----------
    x_extrap = [1, 2, 4] ./ nTrials;

    L = numel(lattice);
    lattice_corr = zeros(1,L);

    % Means across halves/quarters and over shuffles
    mean_half = squeeze(mean(lattice2_shuff, [2 3]));  % L x 1
    mean_quar = squeeze(mean(lattice4_shuff, [2 3]));  % L x 1

    for i = 1:L
        y = [lattice(i), mean_half(i), mean_quar(i)];
        pfit = polyfit(x_extrap, y, 2);   % y = a x^2 + b x + c
        lattice_corr(i) = pfit(3);        % intercept at x = 0
    end
end

% ================= helpers =================

function lat = lat_from_rows_gauss(obs, rows, pid)
    S = cov(obs(rows,:), 1);
    S = (S + S.')/2;
    lat = pid.calculate_latvals(S);
end

function lat = lat_from_rows_disc(obs, rows, nb_vec, pid)
    cnt = accumarray(obs(rows,:), 1, nb_vec);
    p = cnt / sum(cnt(:));
    lat = pid.calculate_latvals(p);
end
