//Kiyotaki and Moore (1997)-Spillovers model with the basic model
//En el presente codigo se utiliza el modelo completo bàsico de Kiyotaki y Moore(1997)
//(secciòn 2) para incorporarle dos sectores productivos de los farmers

//Variables endogena

var k1 k2 k q b1 b2 p1 p2 u y1 y2;

//shock de productividad

varexo epsilon;

parameters  a R v e c kss1 kss2 qss bss1 bss2 pss1 pss2 uss yss1 yss2;

a=0.8;
c=0.2;
R=1.01;
v=2;
e=0.2;

// Estado Estacionario

qss=(R/(R-1))*a;
kss1= 0.5*( qss * (R-1)/R + v);
kss2= 0.5*( qss * (R-1)/R + v);
bss1= kss1*a/(R-1);
bss2= kss2*a/(R-1);
kss=kss1+kss2;
uss = (R-1)*qss/R;
pss1= 1;
pss2= 1;
yss1= (a+c)*kss1;
yss2= (a+c)*kss2;


//EL MODELO

model;

// precio del activo tierra
q = q(+1)/R + (k1 +k2 - v);

// demanda de capital para la producciòn del producto 1
k1 = (1/(q-q(+1)/R))*((a*p1*(1-epsilon)+ q)*k1(-1) - R*b1(-1));

// demanda de capital para la producciòn del producto 2
k2 = (1/(q-q(+1)/R))*((a*p2+ q)*k2(-1) - R*b2(-1));

// capital de los farmers
k=k1+k2;

//Restriccion de deuda en el sector 1
b1 = (1/R)*q(1)*k1;

//Restriccion de deuda en el sector 2
b2 = (1/R)*q(1)*k2;

//funciòn de producciòn del sector 1
y1= (a*(1-epsilon)+c)*k1;

//funciòn de producciòn del sector 2
y2= (a+c)*k2;

//precio del producto1
p1=((a*k1)^(-e))*((((a*(1-epsilon)*k1)^(1-e))+(a*k2)^(1-e))^(e/(1-e)));

//precio del producto 2
p2=((a*k2)^(-e))*((((a*k1)^(1-e))+(a*k2)^(1-e))^(e/(1-e)));

//costo de oportunidad del capital
u = q-q(+1)/R;

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
y1=yss1;
y2=yss2;

end;

shocks;
var epsilon =0.01;

end;

steady;

stoch_simul(hp_filter = 1600, order = 1, irf=20);

