function [ux,uy,a11,a12,a21,a22] = tgv2d(x,y)
% Generates two dimensional Taylor-Green Vortex flow
% INPUT:   
%          - 'x' x-coordinate of input
%          - 'y' y-coordinate of input
%          
% OUTPUT:  
%          - 'ux' x-component of velocity
%          - 'uy' y-component of velocity
%          - 'a11' d ux /dx component of velocity gradient
%          - 'a12' d ux /dy component of velocity gradient
%          -  'a21' d uy /dx component of velocity gradient
%           - 'a22' d uy /dy component of velocity gradient          
% Example:  
% 

%%
ux = cos(x).*sin(y);
uy = -sin(x).*cos(y);
a11 = -sin(x).*sin(y);
a12 = cos(x).*cos(y);
a21 = -cos(x).*cos(y);
a22=-a11;