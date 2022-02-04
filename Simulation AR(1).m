
clc;
clear;
close all;

simulation=10000
% Parámetros del AR(1)
c   = 1;
phi = 0.5;
sigma2 = 1;


%% Simulación para T=50
%==========================
% [1] Simulación
%==========================
% Horizonte temporal
T = 50;
% Generación de secuencia GWN
rng(123); %Semilla
u = randn(T,simulation);  % GWN(0,1)
e = sqrt(sigma2)*u;       % GWN(0, sigma²)

% Aplicación del DGP para el proceso AR(1)
y50 = zeros(T,simulation);
y50(1,:) = c + e(1,:); % Supuesto: y(0) = 0, e(0) = 0
for t = 2:T
    y50(t,:) = c + phi*y50(t-1,:) + e(t,:);
end

% Gráfica de y
figure(1); clf;
plot([y50]); shg; grid on;

%====================
%[2] Estimación OLS
%====================

fprintf('2) Estimación OLS:\n\n');

% Estimación
beta_ols50 = zeros(2,simulation);

for i=1:simulation
X = [ ones(T-1,1) y50(1:T-1,i)];
Y = y50(2:T,i);
beta_ols50(:,i)= ((X'*X)^-1)*X'*Y;
end

% Impresión de la media de los resultados OLS
alpha_hat50= mean(beta_ols50(1, :))';
beta_hat50 = mean(beta_ols50(2,:))';

Tabla_hat = table(alpha_hat50, beta_hat50);
disp(Tabla_hat);

figure(2)
hold on
plot(beta_ols50(1, :));
plot(beta_ols50(2, :));
hold off

%=======================
%[3] Se estima el sesgo
%=======================

%Se estima el sesgo

alpha_bias50=beta_ols50(1, :)-c;
beta_bias50=beta_ols50(2,:)-phi;

figure(3)
hold on
plot(alpha_bias50);
plot(beta_bias50);
hold off
legend('alpha bias', 'beta bias');

 % tabla de la media de los sesgos
 
mean_alpha_bias50= mean(alpha_bias50);
mean_beta_bias50 = mean(beta_bias50)';

Tabla_hat = table(mean_alpha_bias50,mean_beta_bias50);
disp(Tabla_hat);


%% SIMULACIÓN PARA T=250

%==========================
% [1] Simulación
%==========================
% Horizonte temporal
T = 250;
% Generación de secuencia GWN
rng(123); %Semilla
u = randn(T,simulation);  % GWN(0,1)
e = sqrt(sigma2)*u;       % GWN(0, sigma²)
% Aplicación del DGP para el proceso AR(1)
y250 = zeros(T,simulation);
y250(1,:) = c + e(1,:); % Supuesto: y(0) = 0, e(0) = 0
for t = 2:T
    y250(t,:) = c + phi*y250(t-1,:) + e(t,:);
end

% Gráfica de y
figure(1); clf;
plot([y250]); shg; grid on;

%====================
%[2] Estimación OLS
%====================

fprintf('2) Estimación OLS:\n\n');

% Estimación
beta_ols250 = zeros(2,simulation);

for i=1:simulation
X = [ ones(T-1,1) y250(1:T-1,i)];
Y = y250(2:T,i);
beta_ols250(:,i)= ((X'*X)^-1)*X'*Y;
end

% Impresión de la media de los resultados OLS
alpha_hat250= mean(beta_ols250(1, :))';
beta_hat250 = mean(beta_ols250(2,:))';

Tabla_hat = table(alpha_hat250, beta_hat250);
disp(Tabla_hat);

figure(2)
hold on
plot(beta_ols250(1, :));
plot(beta_ols250(2, :));
hold off


%=======================
%[3] Se estima el sesgo
%=======================

%Se estima el sesgo

alpha_bias250=beta_ols250(1, :)-c;
beta_bias250=beta_ols250(2,:)-phi;

figure(3)
hold on
plot(alpha_bias250);
plot(beta_bias250);
hold off
legend('alpha bias', 'beta bias');

 % tabla de la media de los sesgos
 
mean_alpha_bias250= mean(alpha_bias250);
mean_beta_bias250 = mean(beta_bias250)';

Tabla_hat = table(mean_alpha_bias250,mean_beta_bias250);
disp(Tabla_hat);

%% SIMULACIÓN PARA T=500

%==========================
% [1] Simulación
%==========================
% Horizonte temporal
T = 500;

% Generación de secuencia GWN
rng(123); %Semilla
u = randn(T,simulation);  % GWN(0,1)
e = sqrt(sigma2)*u;       % GWN(0, sigma²)

% Aplicación del DGP para el proceso AR(1)
y500 = zeros(T,simulation);

y500(1,:) = c + e(1,:); % Supuesto: y(0) = 0, e(0) = 

for t = 2:T
    y500(t,:) = c + phi*y500(t-1,:) + e(t,:);
end

% Gráfica de y
figure(1); clf;
plot([y500]); shg; grid on;

%====================
%[2] Estimación OLS
%====================

fprintf('2) Estimación OLS:\n\n');

% Estimación
beta_ols500 = zeros(2,simulation);

for i=1:simulation
X = [ ones(T-1,1) y500(1:T-1,i)];
Y = y500(2:T,i);
beta_ols500(:,i)= ((X'*X)^-1)*X'*Y;
end

% Impresión de la media de los resultados OLS
alpha_hat500= mean(beta_ols500(1, :))';
beta_hat500 = mean(beta_ols500(2,:))';

Tabla_hat = table(alpha_hat500, beta_hat500);
disp(Tabla_hat);

figure(2)
hold on
plot(beta_ols500(1, :));
plot(beta_ols500(2, :));
hold off


%=======================
%[3] Se estima el sesgo
%=======================

%Se estima el sesgo

alpha_bias500=beta_ols500(1, :)-c;
beta_bias500=beta_ols500(2,:)-phi;

figure(3)
hold on
plot(alpha_bias500);
plot(beta_bias500);
hold off
legend('alpha bias', 'beta bias');

 % tabla de la media de los sesgos
 
mean_alpha_bias500= mean(alpha_bias500);
mean_beta_bias500 = mean(beta_bias500)';

Tabla_hat = table(mean_alpha_bias500,mean_beta_bias500);
disp(Tabla_hat);

%% SIMULACIÓN PARA T=1000

%==========================
% [1] Simulación
%==========================
% Horizonte temporal
T = 1000;

% Generación de secuencia GWN
rng(123); %Semilla
u = randn(T,simulation);  % GWN(0,1)
e = sqrt(sigma2)*u;       % GWN(0, sigma²)

% Aplicación del DGP para el proceso AR(1)
y1000 = zeros(T,simulation);

y1000(1,:) = c + e(1,:); % Supuesto: y(0) = 0, e(0) = 

for t = 2:T
    y1000(t,:) = c + phi*y1000(t-1,:) + e(t,:);
end

% Gráfica de y
figure(1); clf;
plot([y1000]); shg; grid on;

%====================
%[2] Estimación OLS
%====================

fprintf('2) Estimación OLS:\n\n');

% Estimación
beta_ols1000 = zeros(2,simulation);

for i=1:simulation
X = [ ones(T-1,1) y1000(1:T-1,i)];
Y = y1000(2:T,i);
beta_ols1000(:,i)= ((X'*X)^-1)*X'*Y;
end

% Impresión de la media de los resultados OLS
alpha_hat1000= mean(beta_ols1000(1, :))';
beta_hat1000 = mean(beta_ols1000(2,:))';

Tabla_hat = table(alpha_hat1000, beta_hat1000);
disp(Tabla_hat);

figure(2)
hold on
plot(beta_ols1000(1, :));
plot(beta_ols1000(2, :));
hold off


%=======================
%[3] Se estima el sesgo
%=======================

%Se estima el sesgo

alpha_bias1000=beta_ols1000(1, :)-c;
beta_bias1000=beta_ols1000(2,:)-phi;

figure(3)
hold on
plot(alpha_bias1000);
plot(beta_bias1000);
hold off
legend('alpha bias', 'beta bias');

% tabla de la media de los sesgos
 
mean_alpha_bias1000= mean(alpha_bias1000);
mean_beta_bias1000 = mean(beta_bias1000)';

Tabla_hat = table(mean_alpha_bias1000,mean_beta_bias1000);
disp(Tabla_hat);


%% GRÁFICOS CONJUNTOS

%Histogramas de estimadores de beta por tamaño de muestra
figure(4)
subplot(2,2,1);
histogram(beta_ols50(2,:));
title('T=50');
grid minor ;axis tight;  
set(gca,'MinorGridLineStyle','--','MinorGridAlpha',0.05); 
subplot(2,2,2);
histogram(beta_ols250(2,:));
title('T=250');
grid minor ;axis tight;  
set(gca,'MinorGridLineStyle','--','MinorGridAlpha',0.05); 
subplot(2,2,3);
histogram(beta_ols500(2,:));
title('T=500');
grid minor ;axis tight;  
set(gca,'MinorGridLineStyle','--','MinorGridAlpha',0.05); 
subplot(2,2,4)
histogram(beta_ols1000(2,:));
title('T=1000');
grid minor ;axis tight;  
set(gca,'MinorGridLineStyle','--','MinorGridAlpha',0.05); 


%Histogramas de estimadores de alpha por tamaño de muestra
figure(5)
subplot(2,2,1);
histogram(beta_ols50(1,:));
title('T=50');
grid minor ;axis tight;  
set(gca,'MinorGridLineStyle','--','MinorGridAlpha',0.05); 
subplot(2,2,2);
histogram(beta_ols250(1,:));
title('T=250');
grid minor ;axis tight;  
set(gca,'MinorGridLineStyle','--','MinorGridAlpha',0.05); 
subplot(2,2,3);
histogram(beta_ols500(1,:));
title('T=500');
grid minor ;axis tight;  
set(gca,'MinorGridLineStyle','--','MinorGridAlpha',0.05); 
subplot(2,2,4)
histogram(beta_ols1000(1,:));
title('T=1000');
grid minor ;axis tight;  
set(gca,'MinorGridLineStyle','--','MinorGridAlpha',0.05); 



