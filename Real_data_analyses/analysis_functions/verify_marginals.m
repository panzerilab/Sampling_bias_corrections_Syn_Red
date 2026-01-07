function verify_marginals(p1, p2)

if size(p1)~=size(p2)
    msg = 'Distribution of different sizes';
    error('Prob:DifSize', msg);
end

for i=1:length(size(p1,1))
    stim_cond1 = squeeze(p1(i,:,:));
    stim_cond2 = squeeze(p2(i,:,:));
    margcomp1 = sum(stim_cond1,1)-sum(stim_cond2,1);
    margcomp2 = sum(stim_cond1,2)-sum(stim_cond2,2);

    if sum(margcomp1.^2, "all")+sum(margcomp2.^2, "all")>.01
        disp("Diference between marginals is too big")
        return
    end
end

disp('Marginals match')
end