
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
k2 = 300000;    % N/m
k3 = 40000;     % N/m
k4 = k2;        % N/m

c1 = 4000;  % N*s/m
c2 = 4000;  % N*s/m

l1 = 2;     % m
l2 = 1;     % m

gdl = 4;

dt = 0.0025;
t = 0:dt:50;

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

L = [k2, 0, 0, 0;    % Matriz de participación de carga. Pn = L*P(t)
    0, k4, 0, 0;
    0, 0, 0, 0; 
    0, 0, 0, 0];

zeta_n = [0.0005;
    0.001;
    0.0275;
    0.028];

[FI, Wn] = eig(K, M);
[X, wx] = eig(M\K);
X = fliplr(X);

Mn = diag(X'*M*X);
Kn = diag(X'*K*X);

Wn_vec = diag(Wn); % Los wn son en realidad w^2
Wnd = Wn_vec.^0.5 .*(1-zeta_n.^2).^0.5;


% Descripción de la carga
A = 0.04;  
lambda = 0.35;

v = 10*1000/3600;       % 40 km/h en m/s
wl = 2*pi*v/lambda;     % frec de la carga, función de la vel del auto


%%
% Respuesta en vibraciones libres analíticamente

y = zeros(4, len_t);

% Condiciones iniciales en coordenadas modales
y0 = FI'*M*x0;
yd0 = FI'*M*xd0;

for i = 1:gdl
    
   y(i, :) = (y0(i)*cos(Wnd(i).*t) ...
   + ((yd0(i) + Wn_vec(i)^0.5*zeta_n(i)*y0(i))./Wnd(i)) ...
   .*sin(Wnd(i).*t)).*exp(-zeta_n(i).*Wn_vec(i)^0.5.*t);
    
end

% Transformación a coordenadas geométricas
x = FI*y;

for i = 1:gdl
    
    figure(i)
    plot(t, x(i, :))
    
end

%% 
% Respuesta a carga periódica armónica de forma analítica


beta_n = wl./(Wn_vec.^0.5);
yf = zeros(4, len_t);


P0f = FI'*(A*diag(L));           % Amplitud de carga generalizada

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
y_ode_total = zeros(gdl, len_t);

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
% Respuesta a carga impulsiva mediante la Integral de Duhamel

% a = 0.05;
% b = 0.04;
% 
% v = 20*1000/3600;
% 
% % pico = ((a/2)/v)
% 
% % tri = triangularPulse(0, (a/2)/v, 0.8, t);
% 
% % figure(8)
% % plot(t, tri)
% 
% P = k2*A.*sin(wl*t);
% 
% 
% y_duha = zeros(4, len_t);
% 
% for i=1:gdl
%     
%     yc = P.*cos(Wn_vec(i)^0.5*t);
%     ys = P.*sin(Wn_vec(i)^0.5*t);
% 
%     
%     A_d = zeros(1, len_t);
%     B_d = zeros(1, len_t);
% 
%     decr = exp(-zeta_n(i)*Wn_vec(i)^0.5*dt);
%     
%     for j = 2:len_t
%         
%         A_d(j) = A_d(j-1)*decr + ...
%             (dt/(2*Mn(i)*Wn_vec(i)^0.5))*(yc(j-1)*decr + yc(j));
%         B_d(j) = B_d(j-1)*decr + ...
%             (dt/(2*Mn(i)*Wn_vec(i)^0.5))*(ys(j-1)*decr + ys(j));
%         
%     end
%     
%     aa = A_d.*sin(Wnd(i)*t);
%     bb = B_d.*cos(Wnd(i)*t);
%     y_duha(i, :) = aa - bb;
% 
% end
% 
% for i=1:gdl
%     figure(i)
%     plot(t, y_duha(i, :))
% end


%% 
% Respuesta con el ode45 en coordenadas geométricas


phase = (((l1 + l2)/lambda) - floor((l1 + l2)/lambda))*2*pi;

A = 0.04;  
lambda = 0.35;

v = 200*1000/3600;       % 40 km/h en m/s
wl = 2*pi*v/lambda;     % frec de la carga, función de la vel del auto

xg = [A*sin(wl*(t'));
    A*sin(wl*(t') - phase);
    zeros(1, len_t);
    zeros(1, len_t)];

P0 = L*xg;
% figure(1)
% plot(t, P0(1, :))
% figure(2)
% plot(t, P0(2, :))

ZERO = zeros(gdl, gdl);
I = eye(gdl, gdl);
[t, x_ode] = ode45(@(t, x_ode) ...
    vibforz_geo(t, x_ode, M, C, K, ZERO, I, L, wl, phase, A), t, [x0; xd0]);

% x0 = [0.01; 0; -0.2; 0];
% [t, x_ode] = ode45(@(t, x_ode) ...
%     vibforz_geo(t, x_ode, M, C, K, ZERO, I, L, wl, 0, 0), t, [x0; xd0]);

x_ode = x_ode';

for i=1:gdl
    figure(i)
    plot(t, x_ode(i, :))

end


%% 
% Función para el ode45 en coordenadas modales

function yd = vibforz(t, y, z, w, P0m, wL)

    A = [0, 1;
        -w^2, -2*z*w];
    L = [0;
        P0m*sin(wL*t)];
    yd = A*y + L;

end

%%
% Función para el ode45 en coordenadas geométricas
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