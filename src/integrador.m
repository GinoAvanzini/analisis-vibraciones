
% 
% Proyecto Integrador de Mecánica Vibratoria
% Análisis de Vibraciones en vehículo con simetría axial de 4 GDL
%


clc
clear all
close all

%%
% Datos del sistema

m1 = 25;    % Kg
m2 = m1;    % Kg
m3 = 1000;  % Kg
J = 500;    % Kg*m^2

k1 = 80000;     % N/m
k2 = 300000;     % N/m
k3 = 40000;      % N/m
k4 = k2;        % N/m

c1 = 4000;  % N*s/m
c2 = 4000;  % N*s/m

l1 = 2;     % m
l2 = 1;     % m

gdl = 4;

dt = 0.001;
t = 0:dt:50;

% Condiciones iniciales en coordenadas geométricas
x0 = [0;
    0;
    0.07;
    0.01];

xd0 = [0;
    0;
    0;
    0];

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

L = [k2;    % Matriz de participación de carga. Pn = L*P(t)
    k4;
    0; 
    0];

zeta_n = [0.05;
    0.01;
    0.01;
    0.01];

[FI, Wn] = eig(K, M);
[X, wx] = eig(M\K);
X = fliplr(X);

Mn = diag(X'*M*X);
Kn = diag(X'*K*X);


%%
% Respuesta en vibraciones libres analíticamente


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
% Respuesta a carga periódica armónica de forma analítica

A = 0.04;  
lambda = 0.3;

v = 10*1000/3600;       % 40 km/h en m/s
wl = 2*pi*v/lambda;     % frec de la carga, función de la vel del auto

beta_n = wl./(Wn_vec.^0.5);
yf = zeros(4, length(t));


P0f = FI'*(A*L);           % Amplitud de carga generalizada

D = ((1-beta_n.^2).^2 + (2.*zeta_n.*beta_n).^2).^-0.5;
theta = atan(2.*zeta_n.*beta_n./(1-beta_n.^2));


for i=1:gdl
    
    yf(i, :) = (P0f(i)*D(i)/Kn(i)*sin(wl*t-theta(i)));
    
end

xf = FI*yf;
for i = 1:gdl
    
    figure(i)
    plot(t, xf(i, :))
    
end


%%
% Respuesta a carga armónica con ode45


P0_m = P0f./Mn;
y_ode_total = zeros(gdl, length(t));

for i = 1:gdl
    
    [t, y_ode] = ode45(@(t, y_ode) vibforz(t, y_ode, zeta_n(i), Wn_vec(i)^0.5, P0_m(i), wl), t, [y0(i); yd0(i)]);
    
    y_ode_total(i, :) = y_ode(:, 1);

    figure(2)
    plot(t, y_ode_total(i, :)) % y(:, 2) es la derivada de y(:, 1)
    hold on
    
    
%     figure(3)
%     error__ = x(i, :) - y(:, 1)';
%     plot(t, error__)
%     hold on

end

x_ode = FI*y_ode_total;
figure(5)
plot(t, x_ode(1, :)) % y(:, 2) es la derivada de y(:, 1)
hold on





%% 
% Función para implementar el ode45

function yd = vibforz(t, y, z, w, P0m, wL)

    A = [0, 1;
        -w^2, -2*z*w];
    L = [0;
        P0m*sin(wL*t)];
    yd = A*y + L;

end
