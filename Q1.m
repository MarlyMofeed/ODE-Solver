%%%%%%%%%%%%%%%%% Question one Using FD , BD for some points
clc;
clear all;
close all;
syms  U00 U01 U02 U03 U04 U10 U11 U12 U13 U14 U20 U21 U22 U23 U24 U30 U31 U32 U33 U34 U40 U42 U43 U44 U51 UIR52 UIR41  U53 U54 U41339 U0866
alpha=2-sqrt(3);
peta=alpha;
alpha1=alpha;
peta1=peta;
peta2=1-peta1;
alpha2=1-alpha1;
eq1=U11-0.5==U01;%%%%%%%%
eq2=U12-0.5==U02;%%%%%%%%
eq3=U13-0.5==U03;%%%%%%%%%
eq4=U10==U11;%%%%%%%%%%%%%
eq5=U21+U01+U10+U12-4*U11==0;
eq6=U02+U22+U13+U11-4*U12==0;
eq7=U23+U03+U14+U12-4*U13==0;
eq8=U14==U13;%%%%%%%%%
eq9=U20==U21;%%%%%%%%%%
eq10=U31+U11+U22+U20-4*U21==0;
eq11=U32+U12+U23+U21-4*U22==0;
eq12=U33+U13+U24+U22-4*U23==0;
eq13=U24==U23;%%%%%%%%%%%%
eq14=U30==U31;%%%%%%%%%
eq15=((UIR41+U30)/2)==U40;%%%%%%%%%%%%%%
eq16=UIR41+U21+U32+U30-4*U31==0;
eq17=U42+U22+U33+U31-4*U32==0;
eq18=U43+U23+U34+U32-4*U33==0;
eq19=U34==U33;%%%%%%%%%%%%%%
eq20=2*((U41339/(alpha*(1+alpha)))-(UIR41/(alpha))+(U31/(1+alpha)))+(U42+U40-2*UIR41)==0;
eq21=UIR52+U32+U43+UIR41-4*U42==0;
eq22=U53+U33+U44+U42-4*U43==0;
eq23=U44==U43;%%%%%%%%%%%%%%%%
eq24=-(peta2/peta1)*(((peta1-peta2)/(peta1*peta2))*U0866-((peta2/peta1)*UIR52))==U51;
eq25=2*((U0866/(peta*(1+peta)))-(UIR52/(peta))+(U53/(1+peta)))+(20-2*UIR52+U42)==0;%%%%%@U5,2
eq26=20+U43+U54+UIR52-4*U53==0;
eq27=U54==U53;%%%%%%%%%%
eq28=((alpha2/alpha1)*UIR41-((alpha1-alpha2)/(alpha2*alpha1))*U41339-(alpha1/alpha2)*U51)==0;
eq29=((peta2/peta1)*UIR52-((peta1-peta2)/(peta2*peta1))*U0866-(peta1/peta2)*U51)==0;
eq30=(U03+U14-0.5)/2==U04;%%%%%%%%
eq31=(U01+U10-0.5)/2==U00;%%%%%%%%%
[x,y]=equationsToMatrix([eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9,eq10,eq11,eq12,eq13,eq14,eq15,eq16,eq17,eq18,eq19,eq20,eq21,eq22,eq23,eq24,eq25,eq26,eq27,eq28,eq29,eq30,eq31]);
N=linsolve(x,y);
fprintf('%g\n',abs(N))

