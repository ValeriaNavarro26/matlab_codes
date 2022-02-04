clear all 
[data,names]=xlsread('desempleo_US.xls',1);
%% a) COMPUTACION DE MOVING-AVERAGES PARA LA SERIE

%PARA DATOS TRIMESTRALES, ELEGIMOS LAMBDA=4

T=size(data,1);
filtered_series=zeros(240,8);
for j=2:9
lambda=j;
if rem(lambda,2)==0
k=lambda/2;
smooth=[];
for i=k+1:T-k
 smooth=[smooth;1/lambda*[data(i-k,1)/2+sum(data(i-k+1:i+k-1))+data(i+k,1)/2]];
end
else
 k=(lambda-1)/2;
 smooth=[];
for m=k+1:T-k 
smooth=[smooth;1/lambda*[sum(data(m-k:m+k))]];
end
end
filtered_series(1:end-2*k,lambda-1)=smooth;
end

%COMPARAMOS LOS GRAFICOS PARA DIFERENTES VALORES DE LAMBDA (2 al 5)
figure(1)
subplot(2,2,1)
years=[1961:0.25:2018.75];
plot(years,filtered_series(4:end-5,1),'color',[255 128 102]./255,'LineWidth',2);
hold on 
plot(years,data(5:end-4),'color',[0 76 153]./255,'LineWidth',0.5);
grid minor ;axis tight;  
set(gca,'MinorGridLineStyle','--','MinorGridAlpha',0.05); 
hold off
title('Lambda=2')
xlabel('years')
ylabel('Unemployment rate (%)')
legend('Smooth','Original')

subplot(2,2,2)
plot(years,filtered_series(4:end-5,2),'color',[255 128 102]./255,'LineWidth',2);
hold on 
plot(years,data(5:end-4),'color',[0 76 153]./255,'LineWidth',0.5);
grid minor ;axis tight;  
set(gca,'MinorGridLineStyle','--','MinorGridAlpha',0.05); 
hold off
title('Lambda=3')
xlabel('years')
ylabel('Unemployment rate (%)')
legend('Smooth','Original')


subplot(2,2,3)
plot(years,filtered_series(3:end-6,3),'color',[255 128 102]./255,'LineWidth',2);
hold on 
plot(years,data(5:end-4),'color',[0 76 153]./255,'LineWidth',0.5);
grid minor ;axis tight;  
set(gca,'MinorGridLineStyle','--','MinorGridAlpha',0.05); 
hold off
title('Lambda=4')
xlabel('years')
ylabel('Unemployment rate (%)')
legend('Smooth','Original')

subplot(2,2,4)
plot(years,filtered_series(3:end-6,4),'color',[255 128 102]./255,'LineWidth',2);
hold on 
plot(years,data(5:end-4),'color',[0 76 153]./255,'LineWidth',0.5);
grid minor ;axis tight;  
set(gca,'MinorGridLineStyle','--','MinorGridAlpha',0.05); 
hold off
title('Lambda=5')
xlabel('years')
ylabel('Unemployment rate (%)')
legend('Smooth','Original')

%COMPARAMOS LOS GRAFICOS PARA DIFERENTES VALORES DE LAMBDA (6 al 9)
figure(2)
subplot(2,2,1)
plot(years, filtered_series(2:end-7,5),'color',[255 128 102]./255,'LineWidth',2);
hold on 
plot(years,data(5:end-4),'color',[0 76 153]./255,'LineWidth',0.5);
grid minor ;axis tight;  
set(gca,'MinorGridLineStyle','--','MinorGridAlpha',0.05);   
hold off
title('Lambda=6');
xlabel('years')
ylabel('Unemployment rate (%)')
legend('Smooth','Original')

subplot(2,2,2)
plot(years,filtered_series(2:end-7,6),'color',[255 128 102]./255,'LineWidth',2);
hold on 
plot(years,data(5:end-4),'color',[0 76 153]./255,'LineWidth',0.5);
grid minor ;axis tight;  
set(gca,'MinorGridLineStyle','--','MinorGridAlpha',0.05); 
hold off
title('Lambda=7');
xlabel('years')
ylabel('Unemployment rate (%)')
legend('Smooth','Original')

subplot(2,2,3)
plot(years,filtered_series(1:end-8,7),'color',[255 128 102]./255,'LineWidth',2);
hold on 
plot(years,data(5:end-4),'color',[0 76 153]./255,'LineWidth',0.5);
grid minor ;axis tight;  
set(gca,'MinorGridLineStyle','--','MinorGridAlpha',0.05); 
hold off
title('Lambda=8');
xlabel('years')
ylabel('Unemployment rate (%)')
legend('Smooth','Original')

subplot(2,2,4)
plot(years,filtered_series(1:end-8,7),'color',[255 128 102]./255,'LineWidth',2);
hold on
plot(years,data(5:end-4),'color',[0 76 153]./255,'LineWidth',0.5);
grid minor ;axis tight;  
set(gca,'MinorGridLineStyle','--','MinorGridAlpha',0.05);   
hold off
title('Lambda=9');
xlabel('years')
ylabel('Unemployment rate (%)')
legend('Smooth','Original')

%% COMPUTACION DE LA SERIE DESESTACIONALIZADA
raw_data=data(3:T-2,1);
x=raw_data-filtered_series(1:end-4,3);
%separamos las estaciones para poder hallar las series estimadas de los
%residuos
all_seasons=[]
seasons=[]
for t=1:4
    
for p=0:58
    seasons=[seasons;x(4*p+t)];
end
   all_seasons=[all_seasons seasons]
   seasons=[]
end
x_hat=mean(all_seasons,1)';
s_i=x_hat-mean(x_hat);
%calculamos el componente estacional 
s_i_series=[]
for m=1:59
  s_i_series=[s_i_series; s_i];
end

seasonally_adjusted=raw_data-s_i_series;

figure(3)
plot(raw_data,'color',[255 128 102]./255,'LineWidth',2);
hold on 
plot(seasonally_adjusted,'color',[0 76 153]./255,'LineWidth',1);
grid minor ;axis tight;  
set(gca,'MinorGridLineStyle','--','MinorGridAlpha',0.05);   
hold off
legend('Raw data','Seasonally adjusted')


%% COMPUTACION DEL SUAVIZADO EXPONENCIAL
exp_series=[];
for z=0.1:0.25:0.85
alpha=z;
exp_smooth=[seasonally_adjusted(1,1);zeros(235,1)];
for q=2:236
exp_smooth(q)=exp_smooth(q-1)-alpha*(exp_smooth(q-1)-seasonally_adjusted(q-1));
end
exp_series=[exp_series exp_smooth];
end

year=1960:0.25:2018.75;
vnames{1}='(\alpha=0.1)';
vnames{2}='(\alpha=0.35)';
vnames{3}='(\alpha=0.6)';
vnames{4}='(\alpha=0.85)';

figure(4)
for i=1:4
subplot(2,2,i)
plot(year,exp_series(:,i),'color',[0 76 153]./255,'LineWidth',1)
grid minor ;axis tight;  
set(gca,'MinorGridLineStyle','--','MinorGridAlpha',0.05); 
xlabel('years')
ylabel('Unemployment rate (%)')
title(['Exponential smooth  ',vnames{i}]);

end
