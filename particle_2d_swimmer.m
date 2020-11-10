clear all

L=2;
N=30;
eta=0.2 ;
ranm1 = randn(N,2*N+1);
ranv1=randn(1,N);
ranm2 = randn(N,2*N+1);
ranv2=randn(1,N);
vs=0.1;
%
Np=2000;
%
x=gpuArray(rand(1,Np))*2-1;
y=gpuArray(rand(1,Np))*2-1;
phi = rand(1,Np)*pi-pi/2
scatter(x,y)
%%
dt = 0.001;
T=0.5;
Nt=T/dt;
Lambda =0.5;

%%
tic 
for j=1:Nt
    [ux,uy,a11,a12,a21,a22]=expfourier2d(x,y,L,eta,N,ranm1,ranm2,ranv1,ranv2);
%     ux=0.1*ones(1,Np);
%     uy=0.1*ones(1,Np);
%     a11=0.1*ones(1,Np);
%     a12=0.1*ones(1,Np);
%     a21=0.1*ones(1,Np);
    x = x+dt*(ux+vs*cos(phi));
    y = y+dt*(uy+vs*sin(phi));
    x=mod(x+1,2)-1;
    y=mod(y+1,2)-1;
    o12 = 0.5*(a12-a21) ;
    s12 = 0.5*(a12+a21) ;
    s11 = a11;
    phi=phi+dt*(o12 +Lambda*(s12.*cos(2*phi)-s11.*sin(2*phi)));
    quiver(x,y,cos(phi),sin(phi))
    theta=atan((-s11+sqrt(s11.^2+s12.^2))./s12);
hold on
%quiver(x,y,ux./sqrt(ux.^2+uy.^2),uy./sqrt(ux.^2+uy.^2))
quiver(x,y,cos(theta),sin(theta))
pause(0.1)
clf
j
end
toc
%% Distribution
nx =100;
ndphi=100;
P = zeros(nx,ndphi);
xmin=10^-5;
xmax=0.2;
dphimin=10^-5;
dphimax =1;
%%
tic
for j=1:50
    [ux,uy,a11,a12,a21,a22]=expfourier2d(x,y,L,eta,N,ranm1,ranm2,ranv1,ranv2);
%     ux=0.1*ones(1,Np);
%     uy=0.1*ones(1,Np);
%     a11=0.1*ones(1,Np);
%     a12=0.1*ones(1,Np);
%     a21=0.1*ones(1,Np);
    x = x+dt*(ux+vs*cos(phi));
    y = y+dt*(uy+vs*sin(phi));
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

for l=1:Np
    
        ang = mod(abs(phi(l)-phi(l+1:end)),pi);
        dx = sqrt((x(l)-x(l+1:end)).^2 +(y(l)-y(l+1:end)).^2);
        ang1 =min([ang;pi-ang],[],1);
        
            indx=1+floor((nx-1)*log(dx/xmin)/log(xmax/xmin));
            indang=1+floor((ndphi-1)*log(ang1/dphimin)/log(dphimax/dphimin));
            indices = find(indx>nx | indx<1);
            indx(indices) = [];
            indang(indices) =[];
            indices = find(indang>ndphi | indang<1);
            indx(indices) = [];
            indang(indices) =[];
            indy =sub2ind(size(P), indx, indang);
            P(indy)=P(indy)+1;
        
end
end
toc
dxvec = xmin*exp((0:(1/(nx-1)):1)*log(xmax/xmin));
dphivec = dphimin*exp((0:(1/(ndphi-1)):1)*log(dphimax/dphimin));
%%
P = P./kron(dxvec',dphivec)
%%
[xx,yy]= meshgrid(dxvec,dphivec);

surf(xx,yy,P )
set(gca,'xScale','log')
set(gca,'yScale','log')
set(gca,'zScale','log') %check ordering of xvec and phivec
xlabel('\delta x')
ylabel('\delta \phi')

%%
x1=gpuArray(rand(1,2*Np))*2-1;
y1=gpuArray(rand(1,2*Np))*2-1;

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
quiver(x,y,ux./sqrt(ux.^2+uy.^2),uy./sqrt(ux.^2+uy.^2))
axis equal
lambda1 = abs(o121)./sqrt(s111.^2+s121.^2);
ind1 = find (lambda1/Lambda>1);
ind2 = find (lambda1/Lambda<1);
scat1= scatter(x1(ind1),y1(ind1),'filled','m');
scat2= scatter(x1(ind2),y1(ind2),'filled','y');
alpha(scat1,0.2)
alpha(scat2,0.2)

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
%%
