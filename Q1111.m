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
eq31=2*U10+2*U01-4*U00==1;%%%%%%%U00%%
eq32=U02+2*U11+U00-4*U01==1;%%%%%%%U01
eq33=U03+2*U12+U01-4*U02==1;%%%%%%%U02
eq34=U04+2*U13+U02-4*U03==1;%%%%%%%U03
eq30=2*U03+2*U14-4*U04==1;%%%%%%U04%%
eq1=2*U11+U20+U00-4*U10==0;%%%%%%%%  U10
eq9=2*U21+U30+U10+-4*U20==0;%%%%%%%%%%%%% U20%%%%%%%
eq2=2*U31+U40+U20-4*U30==0;%%%%%%%% U30
eq15=2*UIR41+2*U30-4*U40==0;%%%%%%%%%%%%%%U40
eq5=U01+U21+U12+U10-4*U11==0; %%%%%%% U11
eq6=U02+U22+U13+U11-4*U12==0;%%%%%%%%% U12
eq7=U03+U23+U14+U12-4*U13==0;%%%%%%%% U13
eq3=2*U13+U24+U04-4*U14==0;%%%%%%%%% U14
eq10=U31+U11+U22+U20-4*U21==0;%%%%%%%% U21
eq11=U32+U12+U23+U21-4*U22==0;%%%% U22
eq12=U33+U13+U24+U22-4*U23==0;%%%%%%%%U23
eq13=2*U23+U34+U14+-4*U24==0;%%%%%%%%%% U24%%%%%%%%%%%%
eq16=UIR41+U21+U32+U30-4*U31==0;%%%%%%%%%%U31
eq17=U42+U22+U33+U31-4*U32==0;%%%%%%%%%%U32
eq18=U43+U23+U34+U32-4*U33==0;%%%%%%U33
eq19=2*U33+U44+U24+-4*U34==0;%%%%%%%%%%%%%%U34
eq20=2*((U41339/(alpha*(1+alpha)))-(UIR41/(alpha))+(U31/(1+alpha)))+(U42+U40-2*UIR41)==0;%%%%%%%U41
eq21=UIR52+U32+U43+UIR41-4*U42==0;%%%U42
eq22=U53+U33+U44+U42-4*U43==0;%%%%%%U43
eq23=2*U43+U54+U34+-4*U44==0;%%%%%%%%%%%%%%%%U44
eq24=(peta2/peta1)*(((peta1-peta2)/(peta1*peta2))*U0866-((peta2/peta1)*UIR52))==-U51;%%%%%%%%%%U51
eq25=2*((U0866/(peta*(1+peta)))-(UIR52/(peta))+(U53/(1+peta)))+(20-2*UIR52+U42)==0;%%%%%@U5,2
eq26=20+U43+U54+UIR52-4*U53==0;%%%%%%%%U53
eq27=2*U53+20+U44+-4*U54==0;%%%%%%%%%%U54
eq28=(alpha2/alpha1)*UIR41-((alpha1-alpha2)/(alpha2*alpha1))*U41339-(alpha1/alpha2)*U51==0;%%%%%%%%UIRALPHA
eq29=(peta2/peta1)*UIR52-((peta1-peta2)/(peta2*peta1))*U0866-(peta1/peta2)*U51==0;%%%%%%%%%UIRPETA
[x,y]=equationsToMatrix([eq1,eq2,eq3,eq5,eq6,eq7,eq9,eq10,eq11,eq12,eq13,eq15,eq16,eq17,eq18,eq19,eq20,eq21,eq22,eq23,eq24,eq25,eq26,eq27,eq28,eq29,eq30,eq31,eq32,eq33,eq34]);
N=linsolve(x,y);
fprintf('%g\n',abs(N))

