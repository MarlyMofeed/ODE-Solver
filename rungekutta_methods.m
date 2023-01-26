%% a) Solve this differential equation Numerically using Runge-Kutta’s method for x = 0 to 5 with a step size of 0.2. Note that the initial conditions are y(0)=1 and y’(0)=0..
%%c) Compare and comment.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%if we are dealing with each value of x as if it is the X0 with new h greater than the old one with 0.2%%%%%%%%%%%%%%%%%%%%%%%%%%
function rungekutta_method
% calculates the numerical values of second order ODE using RK
% you have to enter equations of z and y double dash at the end of this file
% you have to enter intial conditions at the beginings of file
clc % clears the command window
clear % clears the workspace
x0=0; % xn value
y0=1; % yn value
z0=0; % y'n value
counter=1;
syms y(x)
Dy = diff(y);
ode = diff(y,x,2) == -1001*Dy-1000*y;
cond1 = y(0) ==1;
cond2 = Dy(0) == 0;
Yexact=dsolve(ode,cond1,cond2);
fprintf(' after solving the equation to get Yexact=\n');
disp(Yexact);
for x1=0:0.2:5
yexact(counter)=(1/999)*(1000*exp((-x1))-exp((-1000*x1))); 
h=x1-x0;
k1=ydash(z0);% F( Zo )
m1=zdash(y0,z0); % G( Y o , Zo ) 
k2=ydash(z0+(m1*(h/2)));% F(Zo +(m1×h/2))
m2=zdash(y0+(k1*(h/2)),z0+(m1*h/2));% G(Y o +(k1×h/2),Zo +(m1×h/2))
k3=ydash(z0+(m2*(h/2))); % G(Y o +(k2×h/2), Zo +(m2×h/2))
m3=zdash(y0+(k2*h/2),z0+(m2*(h/2)));% G(Y o +(k2×h/2), Zo +(m2×h/2))
k4=ydash(z0+m3*h); % F(Y o +(k3×h), Zo +(m3×h))
m4=zdash(y0+k3*h,z0+m3*h);% G(Y o +(k3×h), Zo +(m3×h))
h_kavg=(h/6)*(k1+2*k2+2*k3+k4);
h_mavg=(h/6)*(m1+2*m2+2*m3+m4);
y1(counter)=y0+h_kavg;
z1=z0+h_mavg;
    z0=z1;
    x0=x1;
    y0=y1(counter);
fprintf(' Runge kutta method for %d iteration \n',counter-1);
% disp('Runge kutta method for ');
% disp(counter);
% disp('iteration');
K={'K1';'K2';'K3';'K4';'Kavg*h'};
ValuesK=[k1;k2;k3;k4;h_kavg];
M={'M1';'M2';'M3';'M4';'Mavg*h'};
ValuesM=[m1;m2;m3;m4;h_mavg];
ValueOfY=y1(counter);
ValueOfYDoubleDash=z1;
ValueOfExactY=yexact(counter);
t=table(K,ValuesK,M,ValuesM);
disp(t);
% disp('VALUE OF Y=');
% fprintf('%10.6f \n',y1(counter));
% disp('VALUE OF Y'' =');
% fprintf('%10.6f \n',z1);
% disp('THE EXACT SOLUTION Y=');
% fprintf('%10.6f \n',yexact(counter));
%if (yexact(counter)>1)
    
   erroryWith=(abs((y1(counter)-yexact(counter))/yexact(counter)))*100;
    erroryWithout=(abs((y1(counter)-yexact(counter))))*100;
% disp('THE ERROR PERCENTAGE=');
% fprintf('%10.2e %% \n',errory);
counter=counter+1;
PercentageOfErrorWithoutDividingByYexact=erroryWithout;
PercentageOfErrorWithDividingByYexact=erroryWith;
c=table(ValueOfY,ValueOfYDoubleDash,ValueOfExactY,PercentageOfErrorWithoutDividingByYexact,PercentageOfErrorWithDividingByYexact);
disp(c);
disp('-----------------------------------------------------------------------------------------------------------------------------------------------')
end
n=1;
for i=0:0.2:5
    x1(n)=i;
    n=n+1;
end
plot(x1,y1,'r')
hold on
plot(x1,yexact,'b')
function f=ydash(z)
% enter the equation of y' here
f=z;
function g=zdash(y,z)
% enter the equation of z' here
g=-1001*z-1000*y;