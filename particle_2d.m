clear all

L=2;
N=30;
eta=0.2 ;
ranm1 = randn(N,2*N+1);
ranv1=randn(1,N);
ranm2 = randn(N,2*N+1);
ranv2=randn(1,N);
%
Np=5000;
%
x=gpuArray(rand(1,Np))*2-1;
y=gpuArray(rand(1,Np))*2-1;
phi = rand(1,Np)*pi-pi/2
scatter(x,y)
%%
dt = 0.001;
T=0.01;
Nt=T/dt;
Lambda =1;

%%
tic 
for j=1:Nt
    [ux,uy,a11,a12,a21,a22]=expfourier2d(x,y,L,eta,N,ranm1,ranm2,ranv1,ranv2);
%     ux=0.1*ones(1,Np);
%     uy=0.1*ones(1,Np);
%     a11=0.1*ones(1,Np);
%     a12=0.1*ones(1,Np);
%     a21=0.1*ones(1,Np);
    x = x+dt*ux;
    y = y+dt*uy;
    x=mod(x+1,2)-1;
    y=mod(y+1,2)-1;
    o12 = 0.5*(a12-a21) ;
    s12 = 0.5*(a12+a21) ;
    s11 = a11;
    phi=phi+dt*(o12 +Lambda*(s12.*cos(2*phi)-s11.*sin(2*phi)));
%     quiver(x,y,cos(phi),sin(phi))
%     theta=atan((-s11+sqrt(s11.^2+s12.^2))./s12);
% hold on
% quiver(x,y,cos(theta),sin(theta))
% pause(0.1)
% clf
j
end
toc
%%

%%
scatter3(x,y)

%%
figure(1)
quiver(x,y,cos(phi),sin(phi))
theta=atan((-s11+sqrt(s11.^2+s12.^2))./s12);
hold on
quiver(x,y,cos(theta),sin(theta))
%figure(2)
%
%scatter3(x,y,abs(o12)./sqrt(s12.^2+s11.^2))
%zlim([-2 2])
%% display
%%
x1=gpuArray(rand(1,Np))*2-1;
y1=gpuArray(rand(1,Np))*2-1;

[ux1,uy1,a111,a121,a211,a221]=expfourier2d(x1,y1,L,eta,N,ranm1,ranm2,ranv1,ranv2);
o121 = 0.5*(a121-a211) ;
    s121 = 0.5*(a121+a211) ;
    s111 = a111;
[ux,uy,a11,a12,a21,a22]=expfourier2d(x,y,L,eta,N,ranm1,ranm2,ranv1,ranv2);
    o12 = 0.5*(a12-a21) ;
    s12 = 0.5*(a12+a21) ;
    s11 = a11;
figure(1)
subplot(2,1,1)
quiver(x1,y1,ux1./sqrt(ux1.^2+uy1.^2),uy1./sqrt(ux1.^2+uy1.^2))
axis equal
subplot(2,1,2)
quiver(x,y,cos(phi),sin(phi))
    theta=atan((-s11+sqrt(s11.^2+s12.^2))./s12);
hold on
%quiver(x,y,ux./sqrt(ux.^2+uy.^2),uy./sqrt(ux.^2+uy.^2))
quiver(x,y,cos(theta),sin(theta))
axis equal
lambda1 = abs(o121)./sqrt(s111.^2+s121.^2);
ind1 = find (lambda1/Lambda>1);
ind2 = find (lambda1/Lambda<1);
% scat1= scatter(x1(ind1),y1(ind1),'filled','m');
% scat2= scatter(x1(ind2),y1(ind2),'filled','y');
% alpha(scat1,0.2)
% alpha(scat2,0.2)

%%
