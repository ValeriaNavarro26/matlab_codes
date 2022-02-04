function logpdf = LRM_NIPPosterior(bs2,y, x)
    %% logpdf = LRM_NIPPosterior(bs2,y, X)
    % Computes the log-posterior density of the linear regresion model with the 
    % non-informative prior
    % Inputs:
    %    bs2          [double]  : a ((k+1) x 1) vector with the values of 
    %                             beta at the first k rows and sigma^2 at
    %                             the last row
    %    y            [double]  : a (n x 1) vector with the data of the 
    %                            endogenous variable
    %    X            [double]  : a (n x k) matrix with the data of the 
    %                            exogenous variables
    % Outputs:
    %    logpdf      [cell]   : a value proportional to the log-posterior
    %                           density at sm2
    %
    % Written by Alan Ledesma ©
    [beta, sigma, n, k, xx] = LRM_NIPostMoments(y,x,true);
    b_ols   = beta{1};
    s2_ols  = sigma{2}/(sigma{1}-1);
    
    b  = bs2(1:end-1);
    s2 = bs2(end);
    
    prop2loglike = -(n/2)*log(s2) ...
                   - 1/(2*s2)*( (n-k)*s2_ols+((b-b_ols)'*xx)*(b-b_ols) );
               
    logpdf = prop2loglike;
    
    
end