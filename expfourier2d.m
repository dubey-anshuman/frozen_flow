function [ux,uy,a11,a12,a21,a22] = expfourier2d(x,y,L,eta,N,ranm1,ranm2,ranv1,ranv2,acent) %acent is the center element of a
%%

% ranm1 = randn(N,2*N+1);
% ranv1=randn(1,N);
% ranm2 = randn(N,2*N+1);
% ranv2=randn(1,N);
nvec=-N:1:N;
nmat = repmat(nvec,2*N+1,1);
mmat=flipud(nmat');


xmat=2*pi*i*kron(x,nmat)/L;
ymat=2*pi*i*kron(y,mmat)/L;
%%
a = zeros(2*N+1);
a(1:N,:)= (1/sqrt(2))*ranm1;
a(N+2:end,:)=rot90(a(1:N,:),2);
a(N+1,1:N) = (1/sqrt(2))*ranv1;
a(N+1,N+2:end) = fliplr(a(N+1,1:N));
%%
b = zeros(2*N+1);
b(1:N,:)= (1/sqrt(2))*ranm2;
b(N+2:end,:)=-rot90(b(1:N,:),2);
b(N+1,1:N) = (1/sqrt(2))*ranv2;
b(N+1,N+2:end) = -fliplr(b(N+1,1:N));
a(N+1,N+1)=acent;
amat = gpuArray(repmat(a,1,length(x)));
bmat = repmat(b,1,length(y));
nmat = repmat(nmat,1,length(x));
mmat = repmat(mmat,1,length(y));
%%
%f=real(sum((cos(vecy)+i*sin(vecy)).*exp(-2*pi^2*eta^2*nvec.^2/L^2),2)/L);
ux = reshape(-real(sum(sum(reshape((amat+i*bmat).*mmat.*exp(xmat+ymat).*exp(-2*pi^2*eta^2*(nmat.^2+mmat.^2)/L^2),[2*N+1,2*N+1,length(x)])))*2*pi*i/L^3),[1, length(x)]);
uy = reshape(real(sum(sum(reshape((amat+i*bmat).*nmat.*exp(xmat+ymat).*exp(-2*pi^2*eta^2*(nmat.^2+mmat.^2)/L^2),[2*N+1,2*N+1,length(x)])))*2*pi*i/L^3),[1, length(x)]);
a11 = reshape(real(sum(sum(reshape((amat+i*bmat).*nmat.*mmat.*exp(xmat+ymat).*exp(-2*pi^2*eta^2*(nmat.^2+mmat.^2)/L^2),[2*N+1,2*N+1,length(x)])))*4*pi^2/L^4),[1, length(x)]);
a12 = reshape(real(sum(sum(reshape((amat+i*bmat).*mmat.*mmat.*exp(xmat+ymat).*exp(-2*pi^2*eta^2*(nmat.^2+mmat.^2)/L^2),[2*N+1,2*N+1,length(x)])))*4*pi^2/L^4),[1, length(x)]);
a21 = reshape(-real(sum(sum(reshape((amat+i*bmat).*nmat.*nmat.*exp(xmat+ymat).*exp(-2*pi^2*eta^2*(nmat.^2+mmat.^2)/L^2),[2*N+1,2*N+1,length(x)])))*4*pi^2/L^4),[1, length(x)]);
a22 = -a11;
