function MIX = mutualInformationXYZ(pxyz)
    % Calculate marginal probabilities
    px = squeeze(sum(pxyz, [2, 3])); % Marginal probability density ofX 
    pyz = squeeze(sum(pxyz, 1)); % Joint probability density of Y  and Z 
    % Initialize mutual information
    MIX = 0;
    % Loop through each value of x, y, and z
    for i = 1:size(pxyz, 1)
        for j = 1:size(pxyz, 2)
            for k = 1:size(pxyz, 3)
                if pxyz(i,j,k) ~= 0 && px(i) ~= 0 && pyz(j,k) ~= 0
                    MIX = MIX + pxyz(i,j,k) * log2(pxyz(i,j,k) / (px(i) * pyz(j,k)));
                end
            end
        end
    end
end