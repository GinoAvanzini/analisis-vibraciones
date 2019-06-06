
% 
%% Proyecto Integrador de Mecánica Vibratoria
%% Análisis de Vibraciones en vehículo con simetría axial de 4 GDL
%


clc
clear all
close all

%% Datos del sistema
% 

m1 = 25;    % Kg
m2 = m1;    % Kg
m3 = 1000;  % Kg
J = 500;    % Kg*m^2

k1 = 80000;     % N/m
k2 = 300000;    % N/m
k3 = 40000;     % N/m
k4 = k2;        % N/m

c1 = 4000;  % N*s/m
c2 = 4000;  % N*s/m

l1 = 2;     % m
l2 = 1;     % m

gdl = 4;

dt = 0.0001;    % f = 10kHz
t = 0:dt:5;

len_t = length(t);

% Condiciones iniciales en coordenadas geométricas
x0 = [0.01;
    0;
    -0.02;
    0];

xd0 = [0;
    0;
    0;
    0];

%% Matrices
% 

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

L = [k2, 0, 0, 0;    % Matriz de participación de carga. Pn = L*P(t)
    0, k4, 0, 0;
    0, 0, 0, 0; 
    0, 0, 0, 0];

% zeta_n = [0.0005;
%     0.001;
%     0.0275;
%     0.028];

[FI, Wn] = eig(K, M);
[X, wx] = eig(M\K);
X = fliplr(X);

Mn = diag(X'*M*X);
Kn = diag(X'*K*X);

Wn_vec = diag(Wn); % Los wn son en realidad w^2
% Wnd = Wn_vec.^0.5 .*(1-zeta_n.^2).^0.5;


% Descripción de la carga
A = 0.04;
lambda = 0.35;

v = 100*1000/3600;       % 40 km/h en m/s
wl = 2*pi*v/lambda;     % frec de la carga, función de la vel del auto


%% Respuesta con el ode45 en coordenadas geométricas
% 

clc

phase = (((l1 + l2)/lambda) - floor((l1 + l2)/lambda))*2*pi;

xg = [A*sin(wl*(t));
    A*sin(wl*(t) - phase);
    zeros(1, len_t);
    zeros(1, len_t)];

% for i=1:gdl
%     figure(i)
%     plot(t, xg(i, :))
% 
% end

x0 = [0.010; 0.010; 0.01; 0.1];

P0 = L*xg;

ZERO = zeros(gdl, gdl);
I = eye(gdl, gdl);
% 
% [t, x_ode] = ode45(@(t, x_ode) ...
%     vibforz_geo(t, x_ode, M, C, K, ZERO, I, L, wl, phase, A), t, [x0; xd0]);

[t, x2_ode] = ode45(@(t, x2_ode) ...
    vibforz_geo(t, x2_ode, M, C, K, ZERO, I, L, wl, phase, A), t, [x0; xd0]);

x2_ode = x2_ode';

for i=1:gdl
    figure(i)
    plot(t, x2_ode(i, :))

end

% setGlobalcounter(1);
% [t, x_ode] = ode45(@(t, x_ode) ...
%     vibforz_geom(t, x_ode, M, C, K, ZERO, I, L, xg), t, [x0; xd0]);

% x_ode = x_ode';
% 
% for i=1:gdl
%     figure(i)
%     plot(t, x_ode(i, :))
% 
% end

%% Método de diferencia central para respuesta ante carga general




%% Función para el ode45 en coordenadas modales
% 

function yd = vibforz(t, y, z, w, P0m, wL)

    A = [0, 1;
        -w^2, -2*z*w];
    L = [0;
        P0m*sin(wL*t)];
    yd = A*y + L;

end


%% Función para el ode45 en coordenadas geométricas
%

function x_forz = vibforz_geo(t, x, M, C, K, ZERO, I, L, wl, phase, amp)

    A = [ZERO, I;
        -inv(M)*K, -inv(M)*C];
    
    B = [zeros(4, 1);
        -inv(M)*(L*[amp*sin(wl*t);
                    amp*sin(wl*t - phase);
                    0;
                    0])];
    
    x_forz = A*x + B;

end

% function x_forz = vibforz_geom(t, x, M, C, K, ZERO, I, L, xg)
%     
%     if (getGlobalcounter() > 150001)
% %         x_forz = zeros(8, 1);
%         setGlobalcounter(getGlobalcounter()-1);
%         disp("hola")
%     end
%     
%     a = [ZERO, I;
%         -inv(M)*K, -inv(M)*C];
%     
%     
%     B=vertcat(zeros(4, 1), -inv(M)*(L*xg(:, getGlobalcounter())));
%     
%     x_forz = a*x + B;
% 
%     setGlobalcounter(getGlobalcounter()+1);
%         
% end
% 
% function setGlobalcounter(val)
% global counter
% counter = val;
% 
% end
% 
% function r = getGlobalcounter()
% global counter
% r = counter;
% end

