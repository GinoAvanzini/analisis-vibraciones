
% 
% Proyecto Integrador de Mecánica Vibratoria
% Análisis de Vibraciones en vehículo con simetría axial de 4 GDL
%


clc
clear all
close all


% Datos del sistema

m1 = 100;   % Kg
m2 = m1;    % Kg
m3 = 2000;  % Kg
J = 500;   % Kg*m^2

k1 = 80;    % KN/m
k2 = 50;    % KN/m
k3 = 10;    % KN/m
k4 = k3;    % KN/m

c1 = 10;
c2 = 10;

l1 = 2;     % m
l2 = 1;     % m


% Matrices

M = [m1, 0, 0, 0;
    0, m2, 0, 0;
    0, 0, m3, 0;
    0, 0, 0, J];

C = [c1, 0, -c1, -c1*l1;
    0, c2, -c2, c2*l2;
    -c1, -c2, c1+c2, c1*l1 - c2*l2;
    -c1*l1, c2*l2, c1*l1 - c2*l2, c1*l1^2 + c2*l2^2];
 
K = [k1 + k2, 0, -k1, -k1*l1;
    0, k3 + k4, -k3, k3*l2;
    -k1, -k3, k1 + k3, k1*l1 - k3*l2;
    -k1*l1, k3*l2, k1*l1 - k3*l2, k1*l1^2 + k3*l2^2];


[FI, Wn] = eig(K, M);


%