clc;clear;
%{
.............................................................
Mc=Mass of carbon(Mev/c2).
mp=mass of proton(Mev/c2).
E1=Energy of incident carbon in lab frame.
E3,E4=Energy of outgoing carbon, outgoing proton in lab frame;

E1cm,E2cm= Energy of incident carbon and incident proton in CM Frame.
P1cm,P2cm=Momentum of incident carbon and incident proton in CM Frame.
E3cm,E4cm=Energy of outgoing carbon and outgoing proton in CM Frame.
P3cm,P4cm=Momentum of outgoing carbon and outgoing proton in CM Frame.

thetacm=Angle between outgoing carbon(CM frame) with incident carbon(Lab).
theta3=Angle between outgoing carbon(Lab frame) with incident carbon(Lab).
theta4=Angle between outgoing proton(Lab frame) with incident carbon(Lab).
.............................................................
%}
u=931.49;
mp=938.272;
Mc=10.016853*u;
E3=(Mc+270:1:Mc+390);P3=sqrt((E3.^2)-(Mc^2));
E4=(mp+1:1:mp+130);P4=sqrt((E4.^2)-(mp^2));
E1=400+Mc; P1=sqrt((E1^2)-(Mc^2));
vcm=P1/(E1+mp);
gcm=1/sqrt(1-(vcm^2));

Ecm=sqrt((Mc^2)+(mp^2)+2*E1*mp);
E1cm=((Ecm^2)+(Mc^2)-(mp^2))/(2*Ecm);
P1cm=sqrt((E1cm^2)-(Mc^2));
E4cm=((Ecm^2)+(mp^2)-(Mc^2))/(2*Ecm);
P4cm=sqrt((E4cm^2)-(mp^2));
E3cm=Ecm-E4cm;
P3cm=sqrt((E3cm^2)-(Mc^2));

%plot between theta4 and KE4

Num4=E1*E4+mp*(-E1+E4)-mp*mp;
Den4=sqrt((E1^2)-(Mc^2))*sqrt((E4.^2)-mp*mp);
theta4=acosd(Num4./Den4);
KE4=E4-mp;
plot(theta4,KE4)
xlabel('\theta4')
ylabel('KE4(Mev)')

%.........................................

%plot between theta3 and KE3
%{
Num3=E1*E3+mp*(E3-E1)-Mc*Mc;
Den3=sqrt((E1^2)-(Mc^2))*sqrt((E3.^2)-(Mc^2));
theta3=acosd(Num3./Den3);
KE3=E3-Mc;
plot(theta3,KE3)
xlabel('\theta3')
ylabel('KE3(Mev)')
%}
%..........................................

%plot between thetacm and theta4
%{
thetacm=linspace(0.1,180,20);
theta4=atand(sind(thetacm)./(gcm*( -cosd(thetacm)+(vcm*E4cm/P4cm)) ));
plot(thetacm,theta4)
xlabel('\thetacm')
ylabel('\theta4(Proton)')
%}
%.........................................

%plot between thetacm and theta3
%{
thetacm=linspace(0.1,180,20);
theta3=atand(sind(thetacm)./(gcm*( cosd(thetacm)+(vcm*E3cm/P3cm)) ));
plot(thetacm,theta3)
xlabel('\thetacm)')
ylabel('\theta3(Outgoing_Carbon)')
%}
%...............................