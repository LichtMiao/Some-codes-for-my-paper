syms x1 x2;
f=x1-x2+2*x1^2+2*x1*x2+x2^2;
x=[0;0];
e=10^(-20);
[k ender]=steepest(f,x,e)
