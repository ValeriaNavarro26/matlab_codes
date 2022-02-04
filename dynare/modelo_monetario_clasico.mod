%Modelo Monetario Clásico
%Galí (2008), capítulo 2

%1. Variables

var 
y n w_p pi i delta_m r a
;

%2. Choques

varexo 
e_a
e_m
;

%3. Parámetros

parameters 
sigma psi rho alpha phi_pi phi_a eta
;

sigma=1;
psi=1;
rho=0.9;
alpha=0.33;
phi_pi=1.5;
phi_a=0.9;
eta=4;

%4. Modelo

model;
w_p=sigma*y+psi*n;
y=y(+1)-(1/sigma)*(i-pi(+1)-rho);
w_p=log(1-alpha)+a-alpha*n;
y=a+(1-alpha)*n;
i=rho+phi_pi*pi+e_m;
a=phi_a*a(-1)+e_a;
delta_m=(y-y(-1))-eta*(i-i(-1))+pi;
r=i-pi(+1);
end;

%5. Valores iniciales

initval;
y=0;
n=0;
w_p=0;
pi=0;
a=0;
i=0;
delta_m=0;
r=0;
end;

%6. Solución

steady;
resid;
check;

%7. Magnitud del choque

shocks;
var e_a; stderr 1;
var e_m; stderr 1;
end;

%8. Simulación

stoch_simul(order=1, irf_shocks=(e_a)); 

subplot(3,3,1);
plot(y_e_a, 'k', 'LineWidth', 2.5);
title('Producto (y)', 'fontweight', 'bold', 'fontsize', 15);

subplot(3,3,2);
plot(w_p_e_a, 'k', 'LineWidth', 2.5);
title('Salario real (w-p)', 'fontweight', 'bold', 'fontsize', 15);

subplot(3,3,3);
plot(pi_e_a, 'k', 'LineWidth', 2.5);
title('Inflación (\pi)', 'fontweight', 'bold', 'fontsize', 15);

subplot(3,3,4);
plot(i_e_a, 'k', 'LineWidth', 2.5);
title('Tasa interés nominal (i)', 'fontweight', 'bold', 'fontsize', 15);

subplot(3,3,5);
plot(delta_m_e_a, 'k', 'LineWidth', 2.5);
title('Tasa crecimiento demanda dinero (\Delta m)', 'fontweight', 'bold', 'fontsize', 15);

subplot(3,3,6);
plot(r_e_a, 'k', 'LineWidth', 2.5);
title('Tasa interés real (r)', 'fontweight', 'bold', 'fontsize', 15);

subplot(3,3,7);
plot(a_e_a, 'k', 'LineWidth', 2.5);
title('Productividad (a)', 'fontweight', 'bold', 'fontsize', 15);

stoch_simul(order=1, irf_shocks=(e_m));

subplot(2,1,1);
plot(pi_e_m, 'k', 'LineWidth', 2.5);
title('Inflación (\pi)', 'fontweight', 'bold', 'fontsize', 15);

subplot(2,1,2);
plot(delta_m_e_m, 'k', 'LineWidth', 2.5);
title('Tasa crecimiento demanda dinero (\Delta m)', 'fontweight', 'bold', 'fontsize', 15);