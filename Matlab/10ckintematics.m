clc;clear;
E1=100;
E3=(120000:-1:100000);
u=931.49;
mp=938.272;
Mc=10.0168*u;
theta=acosd((E1*E3+mp*(E3-E1)-Mc*Mc)./( sqrt((E1^2)-(Mc^2))*sqrt((E3.^2-Mc^2)) ))
