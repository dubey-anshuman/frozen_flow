function [f] = expfourier(x,L,eta,N)
%%
nvec = -N:1:N;
vecy=2*pi*kron(nvec,x)/L;
f=real(sum((cos(vecy)+i*sin(vecy)).*exp(-2*pi^2*eta^2*nvec.^2/L^2),2)/L);