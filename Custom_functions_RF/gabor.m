clear
close all
clc

a = 2;
k = 4;

x = -5:0.1:5;
y = -5:0.1:5;

[X,Y] = meshgrid(x,y);

model = makegabor(a, k, 0, -1, -1, 2, X, Y);

figure(1)
contourf(X, Y, model, 10)
colormap(turbo)
colorbar
axis equal

%{
figure(1)
subplot(131)
contourf(X, Y, real(makegabor(sigma, beta, omega, nu, X, Y)), 10)
colormap(turbo)
caxis([-1 1])
axis equal

subplot(132)
contourf(X, Y, imag(makegabor(sigma, beta, omega, nu, X, Y)), 10)
colormap(turbo)
caxis([-1 1])
axis equal

subplot(133)
contourf(X, Y, abs(makegabor(sigma, beta, omega, nu, X, Y)), 10)
colormap(turbo)
caxis([-1 1])
axis equal
%}

function [out] = makegabor(a, k, Ax, Bx, Ay, By, X, Y)

    out = exp(-X.^2/a^2 - Y.^2/a^2).*(Ax*cos(2*pi*X/k) + Bx*sin(2*pi*X/k) + Ay*cos(2*pi*Y/k) + By*sin(2*pi*Y/k));

end




