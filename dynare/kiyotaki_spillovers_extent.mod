
//Kiyotaki and Moore (1997)-Spillovers model with invesment
//En el presente codigo se utiliza el modelo completo de Kiyotaki y Moore(1997)
//para incorporarle dos sectores productivos de los farmers
 

//Variables endogena

var k1 k2 k q b1 b2 p1 p2 u ;

//shock de productividad

varexo epsilon;

parameters  a R v e c kss1 phi pai lambda kss2 qss bss1 bss2 pss1 pss2 uss;

a=0.8;
c=0.2;
R=1.01;
v=2;
pai=0.1;
phi=20;
e=0.2;
lambda=.975;

// Estado Estacionario

qss=(R/(R-1))*(pai*a-(1-lambda)*(1-R+pai*R)*phi)/(lambda*pai+(1-lambda)*(1-R+pai*R));
kss1  = 0.5*(qss * (R-1)/R + v);
kss2  = 0.5*(qss * (R-1)/R + v);
kss=kss1+kss2;
bss1  = ((a + lambda*phi - phi)/(R-1)) * kss1 ;
bss2  = ((a + lambda*phi - phi)/(R-1)) * kss2 ;
uss = (R-1)*qss/R;
pss1= ((a*kss1)^(-e))*((((a*kss1)^(1-e))+(a*kss2)^(1-e))^(e/(1-e)));
pss2= ((a*kss2)^(-e))*((((a*kss1)^(1-e))+(a*kss2)^(1-e))^(e/(1-e)));


//EL MODELO

model;

// precio del activo tierra
q = q(+1)/R + (k - v);

// demanda de capital para la producciòn del producto 1
k1 = (1-pai)*lambda*k1(-1)+(pai/(q+phi-q(+1)/R))*( (p1*a*(1+epsilon) +lambda*phi + q)*k1(-1) - R*b1(-1));

// demanda de capital para la producciòn del producto 2
k2 = (1-pai)*lambda*k2(-1)+(pai/(q+phi-q(+1)/R))*( (p2*a +lambda*phi + q)*k2(-1) - R*b2(-1));

//Restriccion de deuda en el sector 1
b1 = R*b1(-1)+q*(k1 - k1(-1))+phi*(k1 -lambda*k1(-1))-a*(1+epsilon)*k1(-1);

//Restriccion de deuda en el sector 1
b2 = R*b2(-1)+q*(k2 - k2(-1))+phi*(k2 -lambda*k2(-1))-a*k2(-1);

// capital de los farmers
k=k1+k2;

//costo de oportunidad del capital
u = q-q(+1)/R;

// precio del producto 1
p1=((a*k1)^(-e))*((((a*k1)^(1-e))+(a*k2)^(1-e))^(e/(1-e)));

//precio del producto2
p2=((a*k2)^(-e))*((((a*k1)^(1-e))+(a*k2)^(1-e))^(e/(1-e)));

end;

initval;
q=qss;
b1=bss1;
b2=bss2;
k1=kss1;
k2=kss2;
k=kss;
p1=pss1;
p2=pss2;
u=uss;
end;

resid;
steady;
check;

shocks;
var epsilon;
periods 1;
values 0.01;
end;

simul(periods=400);

//gràfica de los impulso respuesta

figure
subplot(421),plot(k1(1:40),'k','LineWidth',1.5),title('k1')
subplot(422),plot(k2(1:40),'k','LineWidth',1.5),title('k2')
subplot(423),plot(b1(1:40),'k','LineWidth',1.5),title('b1')
subplot(424),plot(b2(1:40),'k','LineWidth',1.5),title('b2')
subplot(425),plot(q(1:40),'k','LineWidth',1.5),title('q')
subplot(426),plot(u(1:40),'k','LineWidth',1.5),title('u')
subplot(427),plot(p1(1:40),'k','LineWidth',1.5),title('p1')
subplot(4,2,8),plot(p2(1:40),'k','LineWidth',1.5),title('p2')
