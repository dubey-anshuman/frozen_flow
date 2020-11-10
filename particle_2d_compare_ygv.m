clear all

L=2*pi;
%%

Np=20000;
%
x=gpuArray(rand(1,Np))*L-L/2;
y=gpuArray(rand(1,Np))*L-L/2;
phi = rand(1,Np);
scatter(x,y)

%%
% x(1)=-0.45;
% y(1)=0.48;
x(1)=0.1542;
y(1)=-0.07512;

x(2)=x(1) + randn*0.001;
y(2)=y(1) + randn*0.001;
phi(2) =phi(1) + randn*0.001;
%%
% x =[-0.3 -0.27 -0.25 -0.23 -0.2503 0.01 0.05 -0.05];
% y=[0.1619 0.1619 0.1619 0.1619 -0.4933 -0.4296 -0.4296 -0.4296];
%phi = rand(1,Np)*0.01;
dt = 0.005;
T=1;
Nt=T/dt;
Lambda =10;%0.95/sqrt(2);
vs=0;
%%
rx=[];
ry=[];
orec1 =[];

s1rec1 =[];

s2rec1 =[];
ang1=[];
str = [];
tic 
for j=1:Nt
    [ux,uy,a11,a12,a21,a22]=tgv2d(x,y);
%     ux=0.1*ones(1,Np);
%     uy=0.1*ones(1,Np);
%     a11=0.1*ones(1,Np);
%     a12=0.1*ones(1,Np);
%     a21=0.1*ones(1,Np);
    
     
     
    [k1ux,k1uy,k1a11,k1a12,k1a21,k1a22]=tgv2d(x,y);
    k1ux = k1ux+vs*cos(phi);
    k1uy = k1uy+vs*sin(phi);
    o12 = 0.5*(k1a12-k1a21) ;
    s12 = 0.5*(k1a12+k1a21) ;
    s11 = k1a11;
    k1phi=(o12 +Lambda*(s12.*cos(2*phi)-s11.*sin(2*phi)));
    [k2ux,k2uy,k2a11,k2a12,k2a21,k2a22]=tgv2d(x+dt*(k1ux)/2,y+dt*(k1uy)/2);
        k2ux = k2ux+vs*cos(phi+dt*k1phi/2);
    k2uy = k2uy+vs*sin(phi+dt*k1phi/2);
    o12 = 0.5*(k2a12-k2a21) ;
    s12 = 0.5*(k2a12+k2a21) ;
    s11 = k2a11;
    k2phi=(o12 +Lambda*(s12.*cos(2*(phi+dt*k1phi/2))-s11.*sin(2*(phi+dt*k1phi/2))));
    [k3ux,k3uy,k3a11,k3a12,k3a21,k3a22]=tgv2d(x+dt*(k2ux)/2,y+dt*(k2uy)/2);
    k3ux = k3ux+vs*cos(phi+dt*k2phi/2);
    k3uy = k3uy+vs*sin(phi+dt*k2phi/2);
    o12 = 0.5*(k3a12-k3a21) ;
    s12 = 0.5*(k3a12+k3a21) ;
    s11 = k3a11;
    k3phi=(o12 +Lambda*(s12.*cos(2*(phi+dt*k2phi/2))-s11.*sin(2*(phi+dt*k2phi/2))));
    [k4ux,k4uy,k4a11,k4a12,k4a21,k4a22]=tgv2d(x+dt*(k3ux),y+dt*(k3uy));
    k4ux = k4ux+vs*cos(phi+dt*k3phi);
    k4uy = k4uy+vs*sin(phi+dt*k3phi);
    o12 = 0.5*(k4a12-k4a21) ;
    s12 = 0.5*(k4a12+k4a21) ;
    s11 = k4a11;
    k4phi=(o12 +Lambda*(s12.*cos(2*(phi+dt*k3phi))-s11.*sin(2*(phi+dt*k3phi))));
%     k1 = myfun(t,vartemp);
%         k2 = myfun(t,vartemp+dt*k1/2);
%         k3 = myfun(t,vartemp+dt*k2/2);
%         k4 = myfun(t,vartemp+dt*k3);
%         x = vartemp+(1/6)*dt*(k1+2*k2+2*k3+k4);
%    x = x+dt*(ux+vs*cos(phi));
%      y = y+dt*(uy+vs*sin(phi));
     phi=phi+dt*(1/6)*(k1phi+2*k2phi+2*k3phi+k4phi);
    x = x+(1/6)*dt*(k1ux+2*k2ux+2*k3ux+k4ux);
    y = y+(1/6)*dt*(k1uy+2*k2uy+2*k3uy+k4uy);
    x=mod(x+L/2,L)-L/2;
    y=mod(y+L/2,L)-L/2;
    
    
    rx=[rx; x];
    ry=[ry; y];
    %lt1 = [lt1, (o12/Lambda)/sqrt(s11^2+s12^2)];
    ang1=[ang1;phi];
    orec1 =[orec1;o12];
s1rec1 =[s1rec1;s11];
s2rec1 =[s2rec1;s12];
%str=[str;expfourier2dstream(x,y,L,eta,N,ranm1,ranm2,ranv1,ranv2,acent)];
%     quiver(x,y,cos(phi),sin(phi))
%     theta=atan((-s11+sqrt(s11.^2+s12.^2))./s12);
% hold on
% quiver(x,y,cos(theta),sin(theta))
% pause(0.1)
% clf
j
end
toc
%% Distribution
nx =50;
ndphi=49;
P = zeros(nx,ndphi);
xmin=10^-3;
xmax=1;
dphimin=10^-5;
dphimax =1;
%%
tic
for j=1:2
%      [ux,uy,a11,a12,a21,a22]=tgv2d(x,y);
% 
%     
%      
%      
%     [k1ux,k1uy,k1a11,k1a12,k1a21,k1a22]=tgv2d(x,y);
%     k1ux = k1ux+vs*cos(phi);
%     k1uy = k1uy+vs*sin(phi);
%     o12 = 0.5*(k1a12-k1a21) ;
%     s12 = 0.5*(k1a12+k1a21) ;
%     s11 = k1a11;
%     k1phi=(o12 +Lambda*(s12.*cos(2*phi)-s11.*sin(2*phi)));
%     [k2ux,k2uy,k2a11,k2a12,k2a21,k2a22]=tgv2d(x+dt*(k1ux)/2,y+dt*(k1uy)/2);
%         k2ux = k2ux+vs*cos(phi+dt*k1phi/2);
%     k2uy = k2uy+vs*sin(phi+dt*k1phi/2);
%     o12 = 0.5*(k2a12-k2a21) ;
%     s12 = 0.5*(k2a12+k2a21) ;
%     s11 = k2a11;
%     k2phi=(o12 +Lambda*(s12.*cos(2*(phi+dt*k1phi/2))-s11.*sin(2*(phi+dt*k1phi/2))));
%     [k3ux,k3uy,k3a11,k3a12,k3a21,k3a22]=tgv2d(x+dt*(k2ux)/2,y+dt*(k2uy)/2);
%     k3ux = k3ux+vs*cos(phi+dt*k2phi/2);
%     k3uy = k3uy+vs*sin(phi+dt*k2phi/2);
%     o12 = 0.5*(k3a12-k3a21) ;
%     s12 = 0.5*(k3a12+k3a21) ;
%     s11 = k3a11;
%     k3phi=(o12 +Lambda*(s12.*cos(2*(phi+dt*k2phi/2))-s11.*sin(2*(phi+dt*k2phi/2))));
%     [k4ux,k4uy,k4a11,k4a12,k4a21,k4a22]=tgv2d(x+dt*(k3ux),y+dt*(k3uy));
%     k4ux = k4ux+vs*cos(phi+dt*k3phi);
%     k4uy = k4uy+vs*sin(phi+dt*k3phi);
%     o12 = 0.5*(k4a12-k4a21) ;
%     s12 = 0.5*(k4a12+k4a21) ;
%     s11 = k4a11;
%     k4phi=(o12 +Lambda*(s12.*cos(2*(phi+dt*k3phi))-s11.*sin(2*(phi+dt*k3phi))));
% %     k1 = myfun(t,vartemp);
% %         k2 = myfun(t,vartemp+dt*k1/2);
% %         k3 = myfun(t,vartemp+dt*k2/2);
% %         k4 = myfun(t,vartemp+dt*k3);
% %         x = vartemp+(1/6)*dt*(k1+2*k2+2*k3+k4);
% %    x = x+dt*(ux+vs*cos(phi));
% %      y = y+dt*(uy+vs*sin(phi));
%      phi=phi+dt*(1/6)*(k1phi+2*k2phi+2*k3phi+k4phi);
%     x = x+(1/6)*dt*(k1ux+2*k2ux+2*k3ux+k4ux);
%     y = y+(1/6)*dt*(k1uy+2*k2uy+2*k3uy+k4uy);
%     x=mod(x+L/2,L)-L/2;
%     y=mod(y+L/2,L)-L/2;
% %     quiver(x,y,cos(phi),sin(phi))
% %     theta=atan((-s11+sqrt(s11.^2+s12.^2))./s12);
% % hold on
% % quiver(x,y,cos(theta),sin(theta))
% % pause(0.1)
% % clf
% j
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
[xx,yy]= meshgrid(dphivec,dxvec);

surf(xx,yy,P )

set(gca,'xScale','log')
set(gca,'yScale','log')
set(gca,'zScale','log') %check ordering of xvec and phivec
ylabel('\delta x')
xlabel('\delta \phi')
%%
loglog(dphivec,sum(P(18:22,:),1))
hold on
loglog(dphivec,10^4*dphivec.^-1.2) %close to -1.2 for Lambda =10
%% trajectory plot
ind1=1;
ind2=1;
lamb=orec1./sqrt(s1rec1.^2+s2rec1.^2);
figure(1)
subplot(3,1,1)

plot(orec1(:,ind1:ind2))
hold on
subplot(3,1,2)
plot(lamb(:,ind1:ind2))
subplot(3,1,3)
plot(mod(ang1(:,ind1:ind2),200*pi))
%% separation plot
figure(1)
subplot(3,1,1)

plot(rx(:,1)-rx(:,2))
hold on
subplot(3,1,2)
plot(ry(:,1)-ry(:,2))
subplot(3,1,3)
plot(ang1(:,1)-ang1(:,2))
%%
subplot(3,1,2)
plot(lt)
hold on
%%
plot(lt1)
subplot(3,1,3)
plot(ang-ang(1))
hold on
plot(ang1-ang1(1))
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
x1=gpuArray(rand(1,Np))*L-L/2;
y1=gpuArray(rand(1,Np))*L-L/2;

[ux1,uy1,a111,a121,a211,a221]=tgv2d(x1,y1);
o121 = 0.5*(a121-a211) ;
    s121 = 0.5*(a121+a211) ;
    s111 = a111;
[ux,uy,a11,a12,a21,a22]=tgv2d(x,y);
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
x1=gpuArray(rand(1,4000 ))*2-1;
y1=gpuArray(rand(1,4000))*2-1;

[ux1,uy1,a111,a121,a211,a221]=expfourier2d(x1,y1,L,eta,N,ranm1,ranm2,ranv1,ranv2,acent);
o121 = 0.5*(a121-a211) ;
    s121 = 0.5*(a121+a211) ;
    s111 = a111;
[ux,uy,a11,a12,a21,a22]=tgv2d(x,y);
    o12 = 0.5*(a12-a21) ;
    s12 = 0.5*(a12+a21) ;
    s11 = a11;
    %%
    f1 = expfourier2dstream(x1,y1,L,eta,N,ranm1,ranm2,ranv1,ranv2,acent);
    indplus=find(f1>0);
    indneg=find(f1<0);
figure(1)
quiver(x1,y1,ux1./sqrt(ux1.^2+uy1.^2),uy1./sqrt(ux1.^2+uy1.^2))
hold on
plot(rx,ry,'LineWidth',2)
% s1=scatter(x1(indplus),y1(indplus),'filled','m');
% s2=scatter(x1(indneg),y1(indneg),'filled','y');
% alpha(s1,0.2);
% alpha(s2,0.2);
% plot(rx(:,1),ry(:,1))
% plot(rx(:,2),ry(:,2))
% plot(rx(:,3),ry(:,3))
% plot(rx(:,4),ry(:,4))
% plot(rx(:,5),ry(:,5))
%%
scatter3(x1,y1,f1)
%%
figure(1)
subplot(2,1,1)
plot(rx(:,5:6))
subplot(2,1,2)
plot(ry(:,5:6))

%%
plot(mod(ang1,2*pi))
%%
matty=sqrt(s2rec1.^2+s1rec1.^2);
plot(matty(:,5:8))

%%
plot(rx(:,1:3).^2+ry(:,1:3).^2)

%%
[ux,uy,a11,a12,a21,a22]=tgv2d(x,y);
A=[a11(1),a12(1);a21(1),a22(1)]
[vv,dd]=eig((A+A')/2)
dot(vv(:,1),[ux(1) uy(1)])/norm([ux(1) uy(1)])

dot(vv(:,2),[ux(1) uy(1)])/norm([ux(1) uy(1)])

