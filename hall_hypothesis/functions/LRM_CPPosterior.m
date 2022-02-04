function logpdf = LRM_CPPosterior(bs2,y, x, b_prior, Q_prior, nu_prior, eta_prior)
    %% logpdf = LRM_IPPosterior(bs2,y, X, b_prior, V_prior, aa_prior, bb_prior)
    % Computes the log-posterior density of the linear regresion model with the 
    % informative independent prior
    % Inputs:
    %    bs2          [double]  : a ((k+1) x 1) vector with the values of 
    %                             beta at the first k rows and sigma^2 at
    %                             the last row
    %    y            [double]  : a (n x 1) vector with the data of the 
    %                            endogenous variable
    %    X            [double]  : a (n x k) matrix with the data of the 
    %                            exogenous variables
    %    b_prior      [double]  : a (k x 1) vector with the prior mean of
    %                            beta
    %    Q_prior      [double] : a (k x k) matrix with the prior variance of
    %                            beta divided by sigma^2
    %    nu_prior     [double] : a scalar with the first prior DoF of
    %                            sigma^2
    %    eta_prior    [double] : an scalar with the second prior DoF of
    %                            sigma^2
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
    
    % to generate the loglikelihood of the joint distribution of parameteres 
    % loglikelihood of the data + logbeta_p + logsigma_p
    prop2loglike = -(n/2)*log(s2) ...
                   - 1/(2*s2)*( (n-k)*s2_ols+((b-b_ols)'*xx)*(b-b_ols) );
    % Prior probability for beta.  Notes (5.1)               
	logprobbeta  = log(mvnpdf(b,b_prior,Q_prior));
    % Prior probability for sigma. Notes (5.1)
    logprobsigma = -((nu_prior-k)/2)*log(s2) - (eta_prior/2)/s2;
    
    logpdf = prop2loglike+logprobbeta+logprobsigma;
      
end