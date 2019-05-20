
% 
% Proyecto Integrador de Mecánica Vibratoria
% Análisis de Vibraciones en vehículo con simetría axial de 4 GDL
%


clc
clear all
close all

%%
% Datos del sistema

m1 = 100;   % Kg
m2 = m1;    % Kg
m3 = 2000;  % Kg
J = 500;    % Kg*m^2

k1 = 80;    % KN/m
k2 = 50;    % KN/m
k3 = 10;    % KN/m
k4 = k3;    % KN/m

c1 = 10;
c2 = 10;

l1 = 2;     % m
l2 = 1;     % m

gdl = 4;

dt = 0.05;
t = 0:dt:1000;

% Condiciones iniciales en coordenadas geométricas
x0 = [0.2;
    -0.1;
    0;
    1];

xd0 = [0.5;
    0;
    0.2;
    -2];

%%
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

zeta_n = [0.1;
    0.05;
    0.05;
    0.05];

[FI, Wn] = eig(K, M);

%%
% Respuesta en vibraciones libres analíticamente

Mn = diag(M);

% Los wn son en realidad w^2
Wn_vec = diag(Wn);
Wnd = Wn_vec.*(1-zeta_n).^0.5;

y = zeros(4, length(t));

% Condiciones iniciales en coordenadas modales
y0 = FI'*M*x0;
yd0 = FI'*M*xd0;

for i = 1:gdl
    
   y(i, :) = (y0(i)*cos(Wnd(i)^0.5.*t) ...
   + ((yd0(i) + Wn_vec(i)^0.5*zeta_n(i)*y0(i))./Wnd(i)^0.5) ...
   .*sin(Wnd(i)^0.5.*t)).*exp(-zeta_n(i).*Wn_vec(i)^0.5.*t);
    
end

% Transformación a coordenadas geométricas
x = FI*y;

for i = 1:gdl
    
    figure(i)
    plot(t, x(i, :))
    
end

%% 
% 





%