function [beta, sigma, n, k, xx, ixx, xy] = LRM_NIPostMoments(y,x,FlagMarginal)
    %% [beta_pdf, sigma_pdf, n, k, xx, ixx, xy] = LRM_NIPostMoments(y,X,FlagMarginal)
    % Computes posterior moments of the linear regresion model with the 
    % non informative prior.
    % Inputs:
    %    y [double] : a (n x 1) vector with the data of the endogenous
    %                 variable
    %    X [double] : a (n x k) vector with the data of the exogenous
    %                 variables
    %    FlagMarginal [logical] : if true, moments of the posterior marginal 
    %                             density are computed. Otherwise, moments 
    %                             of the posterior conditional density are
    %                             calculated
    % Outputs:
    %    beta_pdf[cell]  : the first (second or third) element of the cell correspond to the 
    %                      first (second or third) parameter of the posterior distribution 
    %    sigma_pdf[cell] : the first (second or third) element of the cell correspond to the 
    %                      first (second or third) parameter of the posterior distribution 
    %    n      [double] : number of observations
    %    k      [double] : number of explanatory variables
    %    xx     [double] : X'X
    %    ixx    [double] : (X'X)^{-1}
    %    xy     [double] : X'y
    %
    % Written by Alan Ledesma ©
    
    [n,k]  = size(x);
    xx     = x'*x;
    ixx    = xx\eye(k);
    xy     = x'*y;
    b_ols  = ixx*xy;
    u      = y-x*b_ols;
    s2_ols = (u'*u)/(n-k);
    % The posterior marginal density for beta is a t-Student while the
    % associated to sigma is an Inverse Gamma (Koop pp. 20-22)
    if FlagMarginal
        beta = {b_ols, s2_ols*ixx,n-k};
        sigma = {(n-k+2)/2,((n-k)/2)*s2_ols};
    else
    % The posterior condicional density for beta is a Normal while the
    % associated to sigma in an Inverse Gamma (Koop pp. 20-22)
        beta = {b_ols, @(s2)(s2*ixx)};
        sigma = {(n-k+2)/2,((n-k)/2)*s2_ols};        
    end
    
end