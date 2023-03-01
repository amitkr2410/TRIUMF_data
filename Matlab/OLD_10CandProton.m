clear;clc;
u=931.49;
%.............................................%
     
    m1=10.016853*u;    %Projectile
    m2=1.00782503223*u;%Rest
    DE=0;
    m3=10.016853*u+DE;
    m4=1.00782503223*u;     %Lighter one
    E1=m1+60;
    N=300;
    theta4=linspace(0,90,N);%Lighter one

%........................................................................
P1=sqrt((E1^2)-(m1^2));Nmax=N;

for(i=1:N)
    
D(i)=4*P1*P1*(cosd(theta4(i)).^2)-4*((m2+E1)^2);
E=-4*(m2+E1)*(m3*m3-m1*m1-m2*m2-m4*m4-2*E1*m2);
F(i)=-4*(m4*m4)*P1*P1*(cosd(theta4(i))^2)-((m3*m3-m1*m1-m2*m2-m4*m4-2*E1*m2)^2);
   
    if(E*E-4*(D(i)*F(i))<0)
    Nmax=i-1;
        if(Nmax==0)
            disp('No values are possible')
        else
        disp('All the values of theta4 are not possible.')
        disp('Maximum value of theta4 is')
        disp(theta4(Nmax))
        end
    break;
    else
    KE4plus(i)=-m4+((-E+sqrt(E*E-4*(D(i)*F(i))))/(2*D(i)));
    KE4minus(i)=-m4+((-E-sqrt(E*E-4*(D(i)*F(i))))/(2*D(i)));
    end
    
end
if(Nmax==N)
disp('All the values of theta4 are possible');
end
subplot(2,2,1);
plot(theta4(1:Nmax),KE4minus(1:Nmax),theta4(1:Nmax),KE4plus(1:Nmax))
xlabel('\theta4')
ylabel('KE4')

%..............................................%
E4=m4+KE4minus;
P4=real(sqrt((E4.^2)-m4*m4));
theta3=atand((P4.*sind(theta4(1:Nmax)))./(P1-(P4.*cosd(theta4(1:Nmax)))));

A=4*(P1^2)*(cosd(theta3(1:Nmax)).^2)-4*((m2+E1)^2);
B=-4*(m2+E1)*((m4*m4)-(m1*m1)-(m2*m2)-(m3^2)-2*E1*m2);
C=-4*m3*m3*P1*P1*(cosd(theta3(1:Nmax)).^2)-((m4*m4-m1*m1-m2*m2-m3*m3-2*E1*m2)^2);
KE3plus=-m3+((-B+sqrt(B*B-4*(A.*C)))./(2*A));
KE3minus=-m3+((-B-sqrt(B*B-4*(A.*C)))./(2*A));
subplot(2,2,2);
plot(theta3,KE3plus,theta3,KE3minus)
xlabel('\theta3')
ylabel('KE3')
%........................................................................

Ecm=sqrt(m1*m1+m2*m2+2*E1*m2);
E4cm=(Ecm*Ecm-m3*m3+m4*m4)/(2*Ecm);
P4cm=sqrt(E4cm*E4cm-m4*m4);
vcm=P1/(E1+m2);
gcm=(E1+m2)/sqrt(m1*m1+m2*m2+2*m2*E1);
G=gcm*gcm*(tand(theta4(1:Nmax)).^2)+1;
H=-2*vcm*(E4cm/P4cm)*gcm*gcm*(tand(theta4(1:Nmax)).^2);
I=gcm*gcm*(tand(theta4(1:Nmax)).^2)*vcm*vcm*((E4cm/P4cm)^2)-1;
thetacm_plus=acosd((-H+sqrt((H.*H)-4*(G.*I)))./(2*G));
thetacm_minus=acosd((-H-sqrt((H.*H)-4*(G.*I)))./(2*G));
subplot(2,2,3)
plot(theta4(1:Nmax),thetacm_plus,theta4(1:Nmax),thetacm_minus)
xlabel('\theta4')
ylabel('\thetacm')
%........................................................................

E3=m3+KE3plus;
P3=sqrt(((E3).^2)-m3*m3);
Q=m1+m2-sqrt(m4*m4-m1*m1-m2*m2+2*E3*(m2+E1)-2*m2*E1-2*P1*(P3.*cosd(theta3)))-m4;
subplot(2,2,4)
plot(theta3,Q)
xlabel('\theta3')
ylabel('Q_value')
