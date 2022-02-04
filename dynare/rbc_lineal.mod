% Modelo de Ciclos Económicos Reales

%----------------------------------------------------------------
% 0. Cierra todas las ventanas
%----------------------------------------------------------------

close all;

%----------------------------------------------------------------
% 1. Definición de variables y parámetros
%----------------------------------------------------------------

var YY CC KK II HH WW RR a g;
varexo e_a e_g;

parameters

alpha beta delta nu sigma rho gamma phi
y_k c_y i_y g_y
c_k h_k k_y w c y h k
;

%----------------------------------------------------------------
% 2. Calibración
%----------------------------------------------------------------

alpha   =   0.33;
beta    =   0.99;
delta   =   0.023;
nu      =   1.75;
sigma   =   2;
rho     =   0.90;  
gamma  =   0.9;
g_y     =   0.2;
k_y     =   beta*alpha/(1-beta*(1-delta));
y_k     =   k_y^(-1);
c_k     =   y_k-delta-g_y*y_k;
h_k     =   (k_y)^(-1/(1-alpha));
w       =   (1-alpha)*(k_y)^(alpha/(1-alpha));
k       =   (w/((h_k)^nu*(c_k)^sigma))^(1/(sigma+nu));
c       =   c_k*k;
y       =   (k_y)^(-1)*k;
h       =   h_k*k;
c_y     =   c/y;
i_y     =   delta*k/y;
phi   =   alpha*y_k/(alpha*y_k+(1-delta));

%----------------------------------------------------------------
% 3. Modelo
%----------------------------------------------------------------

model; 
YY      =   a +alpha*KK(-1)+(1-alpha)*HH;
YY      =   c_y*CC+i_y*II+g_y*g;
KK      =   delta*II+(1-delta)*KK(-1);
nu*HH   =   WW-sigma*CC;
WW      =   YY-HH;
RR      =   sigma*(CC(+1)-CC)-ln(beta);
RR      =   phi*(YY(+1)-KK)-ln(beta);
a       =   rho*a(-1)+e_a;
g       =   gamma*g(-1)+e_g;
end;

%----------------------------------------------------------------
% 4. Resolución
%----------------------------------------------------------------

initval;
YY      =   0;
CC      =   0;
KK      =   0;
II      =   0; 
HH      =   0;
WW      =   0;
RR      =   0;
a       =   0; 
g       =   0;
end;

shocks;
var e_a =   0.01;
var e_g =   0.01;
end;

steady;

stoch_simul(hp_filter = 1600, order = 1);

