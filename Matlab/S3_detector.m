clc;
clear;
%..................S3_detector_angle_calculation...........................%
D=30; % All distance in cm
r(1:25)=0;
thetamin(1:24)=0;
thetamax(1:24)=0;
thetaavg(1:24)=0;
r(1)=1.1;
r(25)=3.5;
d=(r(25)-r(1))/24;
for(i=1:24)
    r(i)=r(1)+d*(i-1);
    thetamin(i)=atand(r(i)/D);
    thetamax(i)=atand((r(i)+d)/D);
    thetaavg(i)=(thetamin(i)+thetamax(i))/2;
end

%.............................................%
    u=931.494; 
    m1=10.016853*u;     %Projectile Carbon-10
    m2=106.905091*u;    %Rest (107Ag)
    DE1=0;
    DE2=0;
    m3=106.905091*u+DE1;
    m4=10.016853*u+DE2;     % Recoile Excited state
    E1=m1+51.925;
    N=length(thetaavg);
    theta4=thetaavg;
    %N=500;
    %theta4=linspace(0,180,N);%Lighter one

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

%.........Now check among KE4plus and KE4minus which one is allowed......%

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

