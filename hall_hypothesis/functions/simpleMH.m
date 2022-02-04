function [xsim,acceptancerate] = simpleMH(x0,n,b,LTOmega,scale,logpostfun)
    xsim = nan(size(x0,1),n+b);
    aceptance = 1;
    kp1 = size(x0,1);
    for i_ = 1:(n+b)
        % Step 1: Propose starting values for parameters
        if i_==1
            xsim(:,i_) = x0;
        else
            candidate = xsim(:,i_-1);
            
        % Step 2: Simulate candidate multivariate
            candidate = candidate + scale*(LTOmega*randn(kp1,1));
           
        % Step 3: Calculate the candidate likelihood of acceptance
            fx_old = logpostfun(xsim(:,i_-1));           
            fx_new = logpostfun(candidate);
            alpha  = min( fx_new-fx_old,0 );
            
        % Step 4: Simulate u and determine acceptance
            eu = rand;     
            if eu < 1e-20
                u = -1e10;
            else
                u = log(eu);
            end
         
            if u <= alpha
                xsim(:,i_) = candidate;
                aceptance = aceptance+1;
            else
                xsim(:,i_) = xsim(:,i_-1);
                % Note that in the next value of i_ xsim is constant
            end
        end
    end
          % Step 5: r=R+B and drop the first B simulations
    xsim = xsim(:,b+1:(n+b));
    acceptancerate = (aceptance/(n+b))*100;
    
end