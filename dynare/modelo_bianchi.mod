% Versión simple del modelo Javier Biachi
%
% Paul Castillo 
% Lima setiembre 2020

%----------------------------------------------------------------
% 0. Cierra todas las ventanas
%----------------------------------------------------------------

close all;

%----------------------------------------------------------------
% 1. Definición de variables y parámetros
%----------------------------------------------------------------

var  PP_N UU BB CC CC_T YY_N YY_T;


varexo eps1 eps2;

parameters

omega sigma beta k_N k_T rho_N rho_T rr CC_T_SS YY_T_SS BB_SS CC_SS PP_N_SS UU_SS YY_N_SS;


%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------
omega=0.31 ;
sigma=2;
beta=0.91;
k_N=0.32 ;
k_T=0.32;
rho_N=0.9;
rho_T=0.901;
rr=0.04;
YY_T_SS=0.5;
YY_N_SS=1;
CC_T_SS=(1-rr*k_T)/(1+k_N*(1-omega)/omega)*YY_T_SS;
BB_SS=(YY_T_SS-CC_T_SS)/rr;
CC_SS=CC_T_SS^omega*YY_N_SS^(1-omega);
PP_N_SS=(1-omega)/omega*(CC_T_SS/YY_N_SS);
UU_SS=CC_SS^(-sigma)-beta*(1+rr)*CC_SS^(-sigma);

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model;

%----------------------------------------------------------------
% EQ1: TIPO DE CAMBIO REAL
%----------------------------------------------------------------

PP_N=(1-omega)/omega*(CC_T/YY_N);

%----------------------------------------------------------------
% EQ2: ECUACIÓN DE EULER
%----------------------------------------------------------------

UU=omega*CC^(1-sigma)/CC_T-omega*beta*(1+rr)*CC(1)^(1-sigma)/CC_T(1);

%----------------------------------------------------------------
% EQ3: ECUACIÓN DE ENDEUDAMIENTO
%----------------------------------------------------------------

BB+k_N*PP_N*YY_N+k_T*YY_T=0;

%----------------------------------------------------------------
% EQ4: CONSUMO AGREGADO
%----------------------------------------------------------------

CC=CC_T^omega*YY_N^(1-omega);

%----------------------------------------------------------------
% EQ5: BALANZA DE PAGOS
%----------------------------------------------------------------

CC_T=YY_T+BB*(1+rr)-BB(1);

%----------------------------------------------------------------
% EQ6: DINAMICA DEL PRODUCTO TRANSABLE
%----------------------------------------------------------------
log(YY_T/YY_T_SS)=rho_T*log(YY_T(-1)/YY_T_SS)-eps1;

%----------------------------------------------------------------
% EQ7: DINAMICA DEL PRODUCTO NO TRANSABLE
%----------------------------------------------------------------
log(YY_N/YY_N_SS)=rho_N*log(YY_N(-1)/YY_N_SS)+eps2;

end;


initval;

PP_N=PP_N_SS;
UU=UU_SS;
BB=BB_SS;
CC=CC_SS;
CC_T=CC_T_SS;
YY_N=YY_N_SS;
YY_T=YY_T_SS;

end;

shocks;
var eps1 =0.01;
var eps2 =0.01;
end;

steady;

stoch_simul(hp_filter = 1600, order = 2)PP_N CC BB;
