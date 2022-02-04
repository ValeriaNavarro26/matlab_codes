clear;
%Directorio
dfolder='C:\Users\diego\OneDrive\Escritorio\Tarea1\Tarea1\Datos\'; 
addpath('C:\Users\diego\OneDrive\Escritorio\Tarea1\Tarea1\functions');

%% DATOS
dfile=strcat(dfolder,'Datos','.xlsx'); % saca los datos
[data0,names]=xlsread(dfile); %lee el excel
data0(:,1)=(log(data0(:,1))-log(lag0(data0(:,1),1)))*100; %log consumo

T0=15; %training sample
names=names(:,:); 
L=1;
X=[];
for j=0:L 
X=[X lag0(data0(T0+1:end,:),j)];
end
X=packr(X); %Eliminar missing values
y=X(:,1);
X=X(:,2:6);
X=[ones(rows(X),1) X(:,[3 1 4 2 5 ])];
k=6;
% mis seis variables exógenas
%% 1. PRIORS 

b_prior=[0 0.9 0 0 0 0]';  


%% 2.PRIORS DISTRIBUTION

% /////////////////Priors and starting values///////////////////////
x0=[];
for j=0:L 
x0=[x0 lag0(data0(1:T0,:),j) ];
end
x0=packr(x0); %Eliminar missing values
y0=x0(:,1);
x0=x0(:,2:6);
x0=[ones(rows(x0),1) x0(:,[3 1 4 2 5 ])];
beta0=x0\y0;
e0=y0-x0*beta0;
sigmax0=(e0'*e0)/T0;
varbeta0=inv(x0'*x0);


DoFprior = 1+k;
Sm2      = 1/sigmax0; %inverse of Sigma2 
nu_prior  = DoFprior;
eta_prior = DoFprior/(Sm2);
Q_prior  = diag(diag(varbeta0));
std_prior = sqrt(Q_prior*sigmax0);% Same as Q_prior
std_sigma_prior = sqrt((eta_prior/2)^2/((((eta_prior-k+2)/2-1)^2)*((eta_prior-k+2)/2-2)));

%% 3.CONJUGATE POSTERIORS
% Marginal distribution
help LRM_CPPostMoments;
fprintf(1,'\n');

FlagMarginal = true;
[beta_ip, sigma_ip] = LRM_CPPostMoments(y,X,b_prior, Q_prior,...
                      nu_prior, eta_prior,FlagMarginal);

b_cp_mean  = beta_ip{1};
b_cp_std   = sqrt(diag((beta_ip{2}*(beta_ip{3}/(beta_ip{3}-2)))));
s2_cp_mean = sigma_ip{2}/(sigma_ip{1}-1);
s2_cp_std  = sqrt(sigma_ip{2}^2/((sigma_ip{1}-1)^2*(sigma_ip{1}-2)));

for i_ = 1:length(b_cp_mean)
    fprintf(1,'b_cp_%d = %4.2f (%2.2f)\n',i_-1,b_cp_mean(i_),b_cp_std(i_));
end
fprintf(1,'s2_ip  = %4.2f (%2.2f)\n',s2_cp_mean,s2_cp_std);
fprintf(1,'\n');

%% 3.1 GIBBS SAMPLING
% Conditional distribution
FlagMarginal = false;
[beta_nip_cond, sigma_nip_cond] = LRM_CPPostMoments(y,X,b_prior, Q_prior,...
                      nu_prior, eta_prior,FlagMarginal);
% Initilize
R = 10000;
B = ceil(R/10);
betas = nan(k,R+B);
sigma2s = nan(1,R+B);
beta0 = b_prior;
sig20 = 1/Sm2;
betas(:,1) = beta0(:);
sigma2s(1,1) = sig20;
varname = {'beta_1 ','beta_2 ','beta_3 ','beta_4 ','beta_5 ','beta_6 ','sigma^2'};
for r_ = 2:(R+B)
    sig2_cond      = sigma2s(:,r_-1);
    beta_mean_post = beta_nip_cond{1};
    beta_var_post  = beta_nip_cond{2}(sig2_cond);
    betas(:,r_)    = mvnrnd(beta_mean_post,beta_var_post);
    sigma2s(:,r_) = 1/gamrnd(sigma_nip_cond{1},1/sigma_nip_cond{2});
end

% Burning
betas = betas(:,B+1:B+R);
sigma2s = sigma2s(:,B+1:B+R);

% Media y desviación estándar
b_nipGS_mean  = mean(betas,2);
s2_nipGS_mean = mean(sigma2s,2);
b_nipGS_std   = std(betas,1,2);
s2_nipGS_std  = std(sigma2s,1,2);

%Intervalo de confianza 95%
b_GS=[];
for i=1:size(betas,1)
b_GS=[b_GS (prctile(betas(i,:),[97.5 2.5])')];
end

fprintf(1,'\n');
fprintf(2,'Coeff.   \t Prior    \t Conjugate   \n');
for i_ = 1:length(b_cp_mean)   
    fprintf(1,'%s    \t %5.2f    \t %5.2f         \n',...
        varname{i_},b_prior(i_),b_nipGS_mean(i_));    
    fprintf(1,'       \t(%5.2f)   \t(%5.2f)    \n',...
        std_prior(i_),b_nipGS_std(i_));        
end
fprintf(1,'%s       \t %5.2f     \t %5.2f     \n',...
    varname{i_+1},sqrt(1/Sm2),s2_nipGS_mean);    
fprintf(1,'         \t(%5.2f)    \t(%5.2f)     \n',...
    std_sigma_prior,s2_nipGS_std);
fprintf(1,'\n');


%% 3.2Metropolis hasting
% Estimación prior no informativa
FlagMarginal = true;
[beta_nip, sigma_nip, n, k, xx, ixx, xy] = LRM_NIPostMoments(y,X,FlagMarginal);
b_nip_mean  = beta_nip{1};
b_nip_std   = sqrt(diag((beta_nip{2}*(beta_nip{3}/(beta_nip{3}-2)))));
s2_nip_mean = sigma_nip{2}/(sigma_nip{1}-1);
s2_nip_std = sqrt(sigma_nip{2}^2/((sigma_nip{1}-1)^2*(sigma_nip{1}-2)));

%%%%%
findmodefun = @(bs2)(-1*LRM_CPPosterior(bs2,y, X, b_prior, Q_prior, nu_prior, eta_prior));
logpostfun  = @(bs2)(+1*LRM_CPPosterior(bs2,y, X, b_prior, Q_prior, nu_prior, eta_prior));

% Find an innitial candidate considering NI as an initial condition
bs2_0   = [b_nip_mean;s2_nip_mean];
opt     = optimset('display','iter');
[bs2_0,fval,exitflag,output,grad,hessian] = fminunc(findmodefun,bs2_0,opt);
LTOmega = chol(hessian\eye(k+1),'lower');
scale   = 0.82; % set such that acceptance rate is around 35%


R = 100000;     % Draws
B = ceil(R/10);% Burns

[xsim,acceptancerate] = simpleMH(bs2_0,R,B,LTOmega,scale,logpostfun);

% Mean and Standard Deviation by columns
b_condmh_mean  = mean(xsim(1:end-1,:),2);
s2_condmh_mean = mean(xsim(end,:),2);
b_condmh_std   = std(xsim(1:end-1,:),1,2);
s2_condmh_std  = std(xsim(end,:),1,2);

%Intervalo de confianza 95%
b_MH=[];
for i=1:6
b_MH=[b_MH (prctile(xsim(i,:),[97.5 2.5])')];
end


fprintf(1,'\n');
fprintf(2,'Coeff.   \t Prior    \t Non Inform. \t Conjugate   \t Conjugate MH\n');
for i_ = 1:length(b_nip_mean)   
    fprintf(1,'%s  \t %5.2f     \t %5.2f       \t %5.2f     \t %5.2f     \n',...
        varname{i_},b_prior(i_),b_nip_mean(i_),b_nipGS_mean(i_),b_condmh_mean(i_));    
    fprintf(1,'         \t(%5.2f)    \t(%5.2f)   \t(%5.2f)    \t(%5.2f)    \n',...
        std_prior(i_),b_nip_std(i_),b_nipGS_std(i_),b_condmh_std(i_));        
end
fprintf(1,'%s  \t %5.2f     \t %5.2f     \t %5.2f    \t %5.2f     \n',...
    varname{i_+1},sqrt(1/Sm2),s2_nip_mean,s2_nipGS_mean,s2_condmh_mean);    
fprintf(1,'         \t(%5.2f)    \t(%5.2f)    \t(%5.2f)        \t(%5.2f)     \n',...
    std_sigma_prior,s2_nip_std,s2_nipGS_std,s2_condmh_std);
fprintf(1,'\n');

%% 

%% 4. CAMBIO DE Q
Q_prior  = diag(diag(varbeta0))/2;

%% 4.1 GIBBS SAMPLING

% Conditional distribution
FlagMarginal = false;
[beta_nip_cond, sigma_nip_cond] = LRM_CPPostMoments(y,X,b_prior, Q_prior,...
                      nu_prior, eta_prior,FlagMarginal);
% Initilize
R = 100000;
B = ceil(R/10);
betas = nan(k,R+B);
sigma2s = nan(1,R+B);
beta0 = b_prior;
sig20 = 1/Sm2;
betas(:,1) = beta0(:);
sigma2s(1,1) = sig20;
varname = {'beta_1 ','beta_2 ','beta_3 ','beta_4 ','beta_5 ','beta_6 ','sigma^2'};
for r_ = 2:(R+B)
    % get beta|sigma2
    sig2_cond      = sigma2s(:,r_-1);
    beta_mean_post = beta_nip_cond{1};
    beta_var_post  = beta_nip_cond{2}(sig2_cond);
    betas(:,r_)    = mvnrnd(beta_mean_post,beta_var_post);
    % get sigma2|beta
    sigma2s(:,r_) = 1/gamrnd(sigma_nip_cond{1},1/sigma_nip_cond{2});
end

% Burning
betas = betas(:,B+1:B+R);
sigma2s = sigma2s(:,B+1:B+R);

b_nipGS_mean  = mean(betas,2);
s2_nipGS_mean = mean(sigma2s,2);
b_nipGS_std   = std(betas,1,2);
s2_nipGS_std  = std(sigma2s,1,2);

%Intervalo de confianza 95%
b_GS=[];
for i=1:size(betas,1)
b_GS=[b_GS (prctile(betas(i,:),[97.5 2.5])')];
end

fprintf(1,'\n');
fprintf(2,'Coeff.   \t Prior    \t Conjugate   \n');
for i_ = 1:length(b_cp_mean)   
    fprintf(1,'%s    \t %5.2f    \t %5.2f         \n',...
        varname{i_},b_prior(i_),b_nipGS_mean(i_));    
    fprintf(1,'       \t(%5.2f)   \t(%5.2f)    \n',...
        std_prior(i_),b_nipGS_std(i_));        
end
fprintf(1,'%s       \t %5.2f     \t %5.2f     \n',...
    varname{i_+1},sqrt(1/Sm2),s2_nipGS_mean);    
fprintf(1,'         \t(%5.2f)    \t(%5.2f)     \n',...
    std_sigma_prior,s2_nipGS_std);
fprintf(1,'\n');


%% 4.2 Metropolis hasting
FlagMarginal = true;
[beta_nip, sigma_nip, n, k, xx, ixx, xy] = LRM_NIPostMoments(y,X,FlagMarginal);
b_nip_mean  = beta_nip{1};
b_nip_std   = sqrt(diag((beta_nip{2}*(beta_nip{3}/(beta_nip{3}-2)))));
s2_nip_mean = sigma_nip{2}/(sigma_nip{1}-1);
s2_nip_std = sqrt(sigma_nip{2}^2/((sigma_nip{1}-1)^2*(sigma_nip{1}-2)));

%%%%%
findmodefun = @(bs2)(-1*LRM_CPPosterior(bs2,y, X, b_prior, Q_prior, nu_prior, eta_prior));
logpostfun  = @(bs2)(+1*LRM_CPPosterior(bs2,y, X, b_prior, Q_prior, nu_prior, eta_prior));

% Find an innitial candidate considering NI as an initial condition
bs2_0   = [b_nip_mean;s2_nip_mean];
opt     = optimset('display','iter');
[bs2_0,fval,exitflag,output,grad,hessian] = fminunc(findmodefun,bs2_0,opt);
LTOmega = chol(hessian\eye(k+1),'lower');
scale   = 0.82; % set such that acceptance rate is around 35%


R = 10000;     % Draws
B = ceil(R/10);% Burns

[xsim,acceptancerate] = simpleMH(bs2_0,R,B,LTOmega,scale,logpostfun);

% Mean and Standard Deviation by columns
b_condmh_mean  = mean(xsim(1:end-1,:),2);
s2_condmh_mean = mean(xsim(end,:),2);
b_condmh_std   = std(xsim(1:end-1,:),1,2);
s2_condmh_std  = std(xsim(end,:),1,2);

%Intervalo de confianza 95%
b_MH=[];
for i=1:6
b_MH=[b_MH (prctile(xsim(i,:),[97.5 2.5])')];
end

fprintf(1,'\n');
fprintf(2,'Coeff.   \t Prior    \t Non Inform. \t Conjugate   \t Conjugate MH\n');
for i_ = 1:length(b_nip_mean)   
    fprintf(1,'%s  \t %5.2f     \t %5.2f       \t %5.2f     \t %5.2f     \n',...
        varname{i_},b_prior(i_),b_nip_mean(i_),b_nipGS_mean(i_),b_condmh_mean(i_));    
    fprintf(1,'         \t(%5.2f)    \t(%5.2f)   \t(%5.2f)    \t(%5.2f)    \n',...
        std_prior(i_),b_nip_std(i_),b_nipGS_std(i_),b_condmh_std(i_));        
end
fprintf(1,'%s  \t %5.2f     \t %5.2f     \t %5.2f    \t %5.2f     \n',...
    varname{i_+1},sqrt(1/Sm2),s2_nip_mean,s2_nipGS_mean,s2_condmh_mean);    
fprintf(1,'         \t(%5.2f)    \t(%5.2f)    \t(%5.2f)        \t(%5.2f)     \n',...
    std_sigma_prior,s2_nip_std,s2_nipGS_std,s2_condmh_std);
fprintf(1,'\n');




