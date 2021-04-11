clc; 
clear;

N=40;
M=20;
xo=log(2.0/5);
xL=log(7.0/5);
tmax=5.0;
x=linspace(xo,xL,N)';
t=linspace(0,tmax,M);
dx=((xL-xo)/(N-1));

dt=(tmax/(M-1));

r=0.065;
sigma=0.8;
c=r/(0.5*sigma^2);
phi=-0.5+(sqrt(2)/2);

f = @(h)(-c*h);

 %%%%% creating the tridiagonal A %%%%%%
% A1 = full(gallery('tridiag', N, 1, -2, 1)) / (dx^(2));
% A2 = (c-1)*full(gallery('tridiag', N, 1, 0, -1))/ (2*dx);
% I = eye(N);
% A1(1,1)=1; A1(1,2)=0; A1(1,3)=0;
% A1(end,end)=1; A1(end,end-1)=0; A1(end,end-2)=0;
% A2(1,1)=0; A2(1,2)=0; A2(1,3)=0;
% A2(end,end)=0; A2(end,end-1)=0; A2(end,end-2)=0;
 
%A = A1 + A2;
 
%%%%%%%% Initial Condition and Boundary Condition%%%%
U=zeros(N,M);
W=zeros(N,M);
U(1,1)=0; %BC
U(end,1)= (7-5*exp(-c*t(1)))/5; %BC
Uc=max(exp(x(2:N-1))-1,0); %IC
U(:,1)=[U(1,1);Uc;U(end,1)];

for q=2:M
    U(1,q)=0; %BC
    U(end,q)=(7-5*exp(-c*t(q)))/5; %BC 
    U(:,q)=[U(1,q);W(2:N-1,q);U(end,q)];
end
%%%%%%% CN %%%%%%%
 %%%%Coefficients of the tridiagonal system
    e=(((c-1)/(4*dx))-(1/(2*dx^2)))*ones(N,1); 
    g=(-((c-1)/(4*dx))-(1/(2*dx^2)))*ones(N,1);
    f=((1/dt)+(c/2)+(1/dx^2))*ones(N,1);
    f(1)=1;
    g(1)=0;
    f(end)=1;
    e(end)=0;
    [l,u]=tridiagLU(e,f,g);


    %%% Create the tridiag matrix
    ee=diag((-((c-1)/(4*dx))+(1/(2*dx^2)))*ones(N-1,1),-1); 
    gg=diag((((c-1)/(4*dx))+(1/(2*dx^2)))*ones(N-1,1),1);
    ff=diag(((1/dt)-(c/2)-(1/dx^2))*ones(N,1));
    
    B=ff+gg+ee;

    %%%%Loop over time steps
    for m=2:M
        d=B*U(:,m-1); 
        d(1)=U(1,m);
        d(end)=U(end,m);
        U(:,m)=tridiagLUsolve(d,e,l,u,U(:,m-1));
    end

figure
[tt,xx]=meshgrid(t,x);
surf(xx,tt,U)
axis([xo xL 0 5 0 1.4]);
xlabel('x','FontSize',12);
ylabel('t','FontSize',12);
zlabel('v','FontSize',12);
title('European option: Crank-Nicolson method');
