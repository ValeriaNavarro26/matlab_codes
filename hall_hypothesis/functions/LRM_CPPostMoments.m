function [beta, sigma] = LRM_CPPostMoments(y, x, b_prior, Q_prior, nu_prior, eta_prior,FlagMarginal)
    %% [beta_pdf, sigma_pdf] = LRM_CPPostMoments(y, X, b_prior, Q_prior, nu_prior, eta_prior,FlagMarginal)
    % Computes posterior moments of the linear regresion model with the 
    % informative conjugate prior.
    % Inputs:
    %    y            [double] : a (n x 1) vector with the data of the 
    %                            endogenous variable
    %    X            [double] : a (n x k) matrix with the data of the 
    %                            exogenous variables
    %    b_prior      [double] : a (k x 1) vector with the prior mean of
    %                            beta
    %    Q_prior      [double] : a (k x k) matrix with the prior variance of
    %                            beta divided by sigma^2
    %    nu_prior     [double] : a scalar with the first prior DoF of
    %                            sigma^2
    %    eta_prior    [double] : an scalar with the second prior DoF of
    %                            sigma^2
    %    FlagMarginal [logical] : if true, moments of the posterior marginal 
    %                             density are computed. Otherwise, momentos 
    %                             of the posterior conditional density are
    %                             calculated
    % Outputs:
    %    beta   [cell]   : the first (second or third) element of the cell correspond to the 
    %                      first (second or third) parameter of the posterior distribution 
    %    sigma  [cell]   : the first (second or third) element of the cell correspond to the 
    %                      first (second or third) parameter of the posterior distribution 
    %
    % Written by Alan Ledesma ©
    
    [beta, sigma, n, k, xx, ixx, ~] = LRM_NIPostMoments(y,x,FlagMarginal);
    b_ols   = beta{1};
    s2_ols  = sigma{2}/(sigma{1}-1);
    
    % Koop pp. 37
    iQ_prior = Q_prior\eye(k);
    Q_post   = (iQ_prior+xx)\eye(k);
    b_post   = Q_post*(iQ_prior*b_prior+xx*b_ols);
    nu_post  = nu_prior+n;
    aux0     = (b_ols-b_prior)'/(Q_prior+ixx);
    aux0     = aux0*(b_ols-b_prior);
    eta_post = eta_prior + (n-k)*s2_ols + aux0;
    % The posterior marginal density for beta is a t-Student while the
    % associated to sigma is an Inverse Gamma (Koop pp. 37)
    if FlagMarginal
        beta = {b_post, (eta_post/(nu_post-k-1))*Q_post,nu_post-k-1};
        sigma = {(nu_post-k+2)/2,eta_post/2};
    else
    % The posterior condicional density for beta is a Normal while the
    % associated to sigma in an Inverse Gamma (Koop pp. 20-22)
        beta = {b_post, @(s2)(s2*Q_post)};
        sigma = {(nu_post-k+2)/2,eta_post/2};
    end
    
end