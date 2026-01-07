function [rateStim1, rateStim2, sharedRate_Stim] = params_simulations(B, alpha, beta, gamma_flat)
%%
%     rateStim1 = B + alpha*[1,1,0,0];
%     rateStim2 = B + alpha*[1,0,1,0];
    rateStim1 = B + alpha*beta*[1,1,0,0] + alpha*(1-beta)*[1,0,1,0];
    rateStim2 = B + alpha*beta*[1,0,1,0] + alpha*(1-beta)*[1,1,0,0];
    sharedRate_Stim = gamma_flat*[1,1,1,1];
end