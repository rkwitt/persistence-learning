function [testStat, thresh, MMDarr] = pl_mmd(inK, m, n, shuff, alpha)
% MMD compute max.-mean discrepancy two-sample test

    % Compute test statistic
    K  = inK(1:m,1:m);
    L  = inK(m+1:end,m+1:end);
    KL = inK(1:m,m+1:end);
    testStat = 1/m^2*sum(sum(K)) +  ...
        1/n^2*sum(sum(L)) - ...
        2/(n*m)*sum(sum(KL));    

    % Bootstrap H0 stat
    MMDarr = zeros(shuff,1);
    for whichSh=1:shuff
        [~,indShuff] = sort(rand(n+m,1));
        KzShuff = inK(indShuff,indShuff);
        K  = KzShuff(1:m,1:m);
        L  = KzShuff(m+1:end,m+1:end);
        KL = KzShuff(1:m,m+1:end);
        MMDarr(whichSh) = 1/m^2*sum(sum(K)) ...
            + 1/n^2*sum(sum(L)) - ...
            2/(n*m)*sum(sum(KL));        
    end
    
    MMDarr = sort(MMDarr);
    thresh = MMDarr(round((1-alpha)*shuff));
end