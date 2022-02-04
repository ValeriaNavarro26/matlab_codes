function [beta, sigma] = LRM_IPPostMoments(y, x, b_prior, V_prior, aa_prior, bb_prior)
    %% [beta_pdf, sigma_pdf] = LRM_IPPostMoments(y, x, b_prior, V_prior, aa_prior, bb_prior)
    % Computes posterior moments of the linear regresion model with the 
    % informative independent prior (*only conditional pdf are computed*)
    % Inputs:
    %    y            [double]  : a (n x 1) vector with the data of the 
    %                            endogenous variable
    %    X            [double]  : a (n x k) matrix with the data of the 
    %                            exogenous variables
    %    b_prior      [double]  : a (k x 1) vector with the prior mean of
    %                            beta
    %    V_prior      [double]  : a (k x k) matrix with the prior variance of
    %                            beta
    %    aa_prior     [double]  : a scalar with the first prior DoF of
    %                            sigma^2
    %    bb_prior    [double]   : an scalar with the second prior DoF of
    %                            sigma^2
    % Outputs:
    %    beta   [cell]   : the first (second or third) element of the cell correspond to the 
    %                      first (second or third) parameter of the posterior distribution 
    %    sigma  [cell]   : the first (second or third) element of the cell correspond to the 
    %                      first (second or third) parameter of the posterior distribution 
    %
    % Written by Alan Ledesma ©
    [beta, sigma, n, k, xx, ~, xy] = LRM_NIPostMoments(y,x,true);
    b_ols   = beta{1};
    s2_ols  = sigma{2}/(sigma{1}-1);
    
    iV_prior = V_prior\eye(k);
    % Mean and Variance of beta conditional s2. Notes (5.7)-(5.9)
    V_post   = @(s2)( (s2*iV_prior+xx)\eye(k) );
    b_post   = @(s2)( V_post(s2)*( s2*iV_prior*b_prior+xy ) );
    % Mean an Variance of s2 conditional beta. Notes (5.10)-(5-12)
    aa_post_half  = (aa_prior - n)/2;
    bb_post_half  = @(bt)( (bb_prior + (n-k)*s2_ols + ( ((bt - b_ols)'*xx)*(bt - b_ols) ) )/2 );
   
    beta = {b_post, V_post};
    sigma = {aa_post_half,bb_post_half};
    
    
end