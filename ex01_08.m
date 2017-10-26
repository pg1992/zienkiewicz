% Date: Sun Oct  8 23:02:23 -03 2017
% Description:
% 	Solution to the exercise 1.8 from the book
% 	Finite Elements & Approximation by Zienkiewicz.
% 	Equation:
% 
% 	\frac{d^2 \phi}{dx^2} + \phi = x,
% 	\phi = 0 at x = 0,
% 	\frac{d \phi}{dx} + \phi = 0 at x = 1.
%
%	Analytic solution:
%
%	\phi(x) = x - 2*sin(x)/(cos(1)+sin(1))

clear all;
close all;
clc;

Lx = 5;

x0 = 0;
xLx = 1;

dx = (xLx - x0)/Lx;

x = (x0:dx:xLx)';

A = toeplitz([-(2-dx^2) 1 zeros(1, Lx - 1)]);
A(end,:) = [zeros(1, Lx - 2) -1 2*dx 1];
f = [dx^2 * x(2:end) ; 0];

phi = A\f;
phi = [0 ; phi(1:end-1)];

x_c = linspace(0, 1, 1000);
phi_c = x_c - 2*sin(x_c)/(cos(1)+sin(1));
plot(x, phi, 'r-', x_c, phi_c, 'k-');
legend('FEM', 'Analytic', 'location', 'northwest');
