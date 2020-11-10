clear all

L=2;
N=16;
eta=0.2;
ranm1 = randn(N,2*N+1);
ranv1=randn(1,N);
ranm2 = randn(N,2*N+1);
ranv2=randn(1,N);
%%
Np=2;
%%
x=rand(1,Np)*0.02;
y=rand(1,Np)*0.02;
phi = rand(1,Np)*pi-pi/2
scatter(x,y)
%%
dt = 0.001;
T=1;
Nt=T/dt;
Lambda =1.2;
r1=[];
r2=[];
tic
for j=1:Nt
    [ux,uy,a11,a12,a21,a22]=expfourier2d(x,y,L,eta,N,ranm1,ranm2,ranv1,ranv2);
    x = x+dt*ux;
    y = y+dt*uy;
    x=mod(x+1,2)-1;
    y=mod(y+1,2)-1;
    o12 = 0.5*(a12-a21) ;
    s12 = 0.5*(a12+a21) ;
    s11 = a11;
    phi=phi+dt*(o12 + Lambda*(s12.*cos(2*phi)-s11.*sin(2*phi)));
    r1=[r1;x];
    r2=[r2;y];
end
toc
%%
scatter(x,y)
%%
scatter(r1(:,1),r2(:,1),'.b')
hold on
scatter(r1(:,2),r2(:,2),'.r')
%%
quiver(x,y,cos(phi),sin(phi))

hold on
%quiver(x,y,cos(theta),sin(theta))

%%
