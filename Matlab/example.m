u=931.49;
mp=938.272;
Mc=10.0168*u;
E1=400+Mc; P1=sqrt((E1^2)-(Mc^2));
E4=mp+20;%(mp+1:1:mp+130);P4=sqrt((E4.^2)-(mp^2));
theta4=linspace(0,1.5,20);
A=E4*(1+(mp/Mc));
C=-(2/Mc)*sqrt(Mc*mp*E1*E4)*cos(theta4);
Q=A+C;
plot(theta4,Q)
ylabel('Q')
xlabel('\theta 4')