clear;clc;
u=931.494; 
%........General..Program....................%
%.............................................%
     
    m1=10.016853*u;     %Projectile
    m2=106.905091*u;    %Rest (107Ag)
    DE1=0;
    DE2=0;
    m3=106.905091*u+DE1;
    m4=10.016853*u+DE2;     % Recoile Excited state
    E1=m1+51.925;
    N=500;
    theta4=linspace(0,180,N);%Lighter one

%........................................................................%
%..............Part-one..Calculate..KE4..of..the..m4...particle...........%
P1=sqrt((E1^2)-(m1^2));Nmax=N;

for(i=1:N)
    
D(i)=4*P1*P1*(cosd(theta4(i)).^2)-4*((m2+E1)^2);
E=-4*(m2+E1)*(m3*m3-m1*m1-m2*m2-m4*m4-2*E1*m2);
F(i)=-4*(m4*m4)*P1*P1*(cosd(theta4(i))^2)-((m3*m3-m1*m1-m2*m2-m4*m4-2*E1*m2)^2);
 
    if(E*E-4*(D(i)*F(i))<0 )
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

%.........Now check among KE4plus and KE4minus, which one is allowed......%

P4_p=sqrt(((m4+KE4plus).^2)-m4*m4);
P4_m=sqrt(((m4+KE4minus).^2)-m4*m4);
E4_p=m4+KE4plus;
E4_m=m4+KE4minus;
P3_p=sqrt(P1*P1+P4_p.^2-2*P1*P4_p.*cosd(theta4(1:Nmax)));
P3_m=sqrt(P1*P1+P4_m.^2-2*P1*P4_m.*cosd(theta4(1:Nmax)));
E3_p=sqrt(P3_p.^2+m3*m3);
E3_m=sqrt(P3_m.^2+m3*m3);
Etotal_p=E3_p+E4_p;
Etotal_m=E3_m+E4_m;
ET=E1+m2;
Erel_p=abs(ET-Etotal_p)*100/ET;
Erel_m=abs(ET-Etotal_m)*100/ET;

BOTH=0;
for(i=1:Nmax)
    if(Erel_m(i)<0.000001 && Erel_p(i)<0.000001 && KE4plus(i)>0&& KE4minus(i)>0 )
        BOTH=1;
        
    else
    if(Erel_p(i)<0.000001 && KE4plus(i)>0)
        KE4(i)=KE4plus(i);
        BOTH=2;
    end
    if(Erel_m(i)<0.000001 && KE4minus(i)>0)
        KE4(i)=KE4minus(i);
        BOTH=3;
    end
    end

         
    
end
%Nmax=length(KE4);
%....................................................................%

subplot(2,2,1);
if(BOTH==1)
plot(theta4(1:Nmax),KE4plus(1:Nmax),'g',theta4(1:Nmax),KE4minus(1:Nmax),'b')
xlabel('\theta4')
ylabel('KE4')
end
if(BOTH==2 || BOTH==3)
    Nmax=length(KE4);
plot(theta4(1:Nmax),KE4(1:Nmax))
xlabel('\theta4')
ylabel('KE4')
end


%.....................................................................%


...............Plot_of_KE3___Vs___Theta3 (If BOTH==1)..............................
............... When we have two good solution for E4 energy............................
if(BOTH==1) 
E4_p=m4+KE4plus;
E4_m=m4+KE4minus;
P4_p=real(sqrt((E4_p.^2)-m4*m4));
P4_m=real(sqrt((E4_m.^2)-m4*m4));
theta3_p=atand((P4_p.*sind(theta4(1:Nmax)))./(P1-(P4_p.*cosd(theta4(1:Nmax)))));
theta3_m=atand((P4_m.*sind(theta4(1:Nmax)))./(P1-(P4_m.*cosd(theta4(1:Nmax)))));
for(i=1:Nmax)
    if((P4_p(i)*sind(theta4(i)))==0 && (P1-(P4_p(i).*cosd(theta4(i))))==0)
    theta3_p(i)=90;  %when theta4=0 theta3 should be 90 degree
    end
    if((P4_m(i)*sind(theta4(i)))==0 && (P1-(P4_m(i).*cosd(theta4(i))))==0)
    theta3_m(i)=90;  %when theta4=0 theta3 should be 90 degree
    end
    if((P4_p(i)*sind(theta4(i)))==0 && (P1-(P4_p(i).*cosd(theta4(i))))>0)
    theta3_p(i)=0;  %
    end
    if((P4_m(i)*sind(theta4(i)))==0 && (P1-(P4_m(i).*cosd(theta4(i))))<0)
    theta3_m(i)=180;  %
    end
    if(theta3_p(i)<0)
        theta3_p(i)=theta3_p(i)+180;
    end
    if(theta3_m(i)<0)
        theta3_m(i)=theta3_m(i)+180;
    end
end
P3_p=sqrt(P1*P1+P4_p.^2-2*P1*P4_p.*cosd(theta4(1:Nmax)));
P3_m=sqrt(P1*P1+P4_m.^2-2*P1*P4_m.*cosd(theta4(1:Nmax)));
E3_p=sqrt(P3_p.^2+m3*m3);
E3_m=sqrt(P3_m.^2+m3*m3);
KE3_p=E3_p-m3;
KE3_m=E3_m-m3;

subplot(2,2,2);
plot(theta3_p,KE3_p,'g',theta3_m,KE3_m,'b')
xlabel('\theta3')
ylabel('KE3')
end
...............Plot_of_KE3___Vs___Theta3 (If BOTH==2 or 3)..............................
...............Only One type of KE4_energy solution is possible.........................
if(BOTH==2 || BOTH==3)
E4=m4+KE4;
P4=real(sqrt((E4.^2)-m4*m4));
theta3=atand((P4.*sind(theta4(1:Nmax)))./(P1-(P4.*cosd(theta4(1:Nmax)))));
for i=1:Nmax
    if((P4(i)*sind(theta4(i)))==0 && (P1-(P4(i)*cosd(theta4(i))))==0)
    theta3(i)=90;  %when theta4=0 theta3 should be 90 degree
    end
    if((P4(i)*sind(theta4(i)))==0 && (P1-(P4(i)*cosd(theta4(i))))>0)
    theta3(i)=0;  %
    end    
    if((P4(i)*sind(theta4(i)))==0 && (P1-(P4(i)*cosd(theta4(i))))<0)
    theta3(i)=180;  %
    end
    if(theta3(i)<0)
        theta3(i)=theta3(i)+180;
    end
end

P3=sqrt(P1*P1+P4.^2-2*P1*P4.*cosd(theta4(1:Nmax)));
E3=sqrt(P3.^2+m3*m3);
KE3=E3-m3;

subplot(2,2,2);
plot(theta3,KE3,'g')
xlabel('\theta3')
ylabel('KE3')
end
%............................Plot..of..Thetacm..vs..theta4.........................%
Ecm=sqrt(m1*m1+m2*m2+2*E1*m2);
E4cm=(Ecm*Ecm-m3*m3+m4*m4)/(2*Ecm);
P4cm=sqrt(E4cm*E4cm-m4*m4);
vcm=P1/(E1+m2);
gcm=(E1+m2)/sqrt(m1*m1+m2*m2+2*m2*E1);
G=gcm*gcm*(tand(theta4(1:Nmax)).^2)+1;
H=-2*vcm*(E4cm/P4cm)*gcm*gcm*(tand(theta4(1:Nmax)).^2);
I=gcm*gcm*(tand(theta4(1:Nmax)).^2)*vcm*vcm*((E4cm/P4cm)^2)-1;
for(i=1:Nmax)
    if(H(i)*H(i)-4*G(i)*I(i)>=0)
        thetacm_p(i)=acosd((-H(i)+sqrt((H(i)*H(i))-4*(G(i)*I(i))))/(2*G(i)));
        thetacm_m(i)=acosd((-H(i)-sqrt((H(i)*H(i))-4*(G(i)*I(i))))/(2*G(i)));
    end
end 

for(i=1:Nmax)
    if(isreal(thetacm_p(i))==1)
    Err_p(i)=abs(tand(theta4(i))-((sind(thetacm_p(i)))./(gcm*(-cosd(thetacm_p(i))+(vcm*E4cm/P4cm)))));
    else
        Err_p(i)=-1;
    end
    if(isreal(thetacm_m(i))==1)
    Err_m(i)=abs(tand(theta4(i))-((sind(thetacm_m(i)))./(gcm*(-cosd(thetacm_m(i))+(vcm*E4cm/P4cm)))));
    else
        Err_m(i)=-1;
    end
end
flag=0;
for(i=1:Nmax)
        if(Erel_p(i)<0.000001 && Erel_m(i)<0.000001)
        flag=1;
        else     
            if(Erel_p(i)<0.000001)
            thetacm(i)=thetacm_p(i);
            flag=2;
            end
            if(Erel_m(i)<0.000001)
            thetacm(i)=thetacm_m(i);
            flag=3;
            end
        end
end
subplot(2,2,3)
if(flag==1)
plot(theta4(1:Nmax),thetacm_p(1:Nmax),'g',theta4(1:Nmax),thetacm_m(1:Nmax),'b')
xlabel('\theta4')
ylabel('\thetacm')
end
if(flag==2 || flag==3)
plot(theta4(1:Nmax),thetacm(1:Nmax),'g')
xlabel('\theta4')
ylabel('\thetacm')
end

%...............................................................%


.................PLOT..of..Q-value..............................

%Q1=m1+m2-sqrt(m4*m4-m1*m1-m2*m2+2*E3*(m2+E1)-2*m2*E1-2*P1*(P3.*cosd(theta3)))-m4;
if(BOTH==1)
Q_p=m1+m2-sqrt(m3*m3-m1*m1-m2*m2+2*E4_p*(m2+E1)-2*m2*E1-2*P1*(P4_p.*cosd(theta4(1:Nmax))))-m3;
Q_m=m1+m2-sqrt(m3*m3-m1*m1-m2*m2+2*E4_m*(m2+E1)-2*m2*E1-2*P1*(P4_m.*cosd(theta4(1:Nmax))))-m3;
subplot(2,2,4)
plot(theta4(1:Nmax),Q_p,'g',theta4(1:Nmax),Q_m,'g')
xlabel('\theta3')
ylabel('Q_value')
end
if(BOTH==2 || BOTH==3)
Q=m1+m2-sqrt(m3*m3-m1*m1-m2*m2+2*E4*(m2+E1)-2*m2*E1-2*P1*(P4.*cosd(theta4(1:Nmax))))-m3;
subplot(2,2,4)
plot(theta4(1:Nmax),Q,'g')
xlabel('\theta3')
ylabel('Q_value')
end
%ylim([-10 +10 ])