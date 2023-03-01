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
E1=400+Mc; P1=sqrt((E1^2)-(Mc^2));

%{
E4=mp+20;P4=sqrt((E4^2)-(mp^2));
theta4=(0:1:30);
Mcn=sqrt((Mc^2)+2*(mp^2)+ 2*E1*mp -2*E4*(mp+E1)+2*P1*P4*cos(theta4) );
Qvalue=(Mc - Mcn);
plot(theta4,Qvalue)
ylabel('Qvalue')
xlabel('\theta4(proton)')
%}


E4=(mp+1:1:mp+130);P4=sqrt((E4.^2)-(mp^2));
theta4=10;
Mcn=sqrt((Mc^2)+2*(mp^2)+ 2*E1*mp -2*E4*(mp+E1)+2*P1*P4*cosd(theta4) );
Qvalue=(Mc - Mcn);
KE4=E4-mp;
plot(KE4,Qvalue)
ylabel('Qvalue')
xlabel('KE4(proton)')
