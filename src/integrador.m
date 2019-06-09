
% 
%% Proyecto Integrador de Mecánica Vibratoria

%% Análisis de Vibraciones en vehículo con simetría axial de 4 GDL
%


clc
clear variables
close all

%% Datos del sistema
% 

m1 = 50;    % Kg
m2 = m1;    % Kg
m3 = 1000;  % Kg
J = 500;    % Kg*m^2

k1 = 80000;     % N/m
k2 = 300000;    % N/m
k3 = 40000;     % N/m
k4 = k2;        % N/m

c1 = 4000;  % N*s/m
c2 = 4000;  % N*s/m

l1 = 1.4;     % m
l2 = 2;     % m

gdl = 4;

dt = 0.0001;    % f = 10kHz
t = 0:dt:5;

len_t = length(t);

% Condiciones iniciales en coordenadas geométricas
x0 = [0
    0;
    0;
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




%% Respuesta con el ode45 en coordenadas geométricas en vibraciones libres
% 

x0 = [0; 0; 0.05; 0];
ZERO = zeros(gdl, gdl);
I = eye(gdl, gdl);

[t, x_ode] = ode45(@(t, x2_ode) ...
    vibforz_geo(t, x2_ode, M, C, K, ZERO, I, L, 0, 0, 0), t, [x0; xd0]);

x_ode = x_ode';

for i=1:gdl
    figure(i)
    plot(t, x_ode(i, :))
    hold on

end

clear x_ode
clear xg

close all

%% Respuesta con el ode45 en coordenadas geométricas en vibraciones 
% forzadas con carga armónica
%

clc

% Descripción de la carga
A = 0.04;
lambda = 0.35;

v = 100*1000/3600;       % 40 km/h en m/s
wl = 2*pi*v/lambda;     % frec de la carga, función de la vel del auto

phase = (((l1 + l2)/lambda) - floor((l1 + l2)/lambda))*2*pi;

% xg = [A*sin(wl*(t));
%     A*sin(wl*(t) - phase);
%     zeros(1, len_t);
%     zeros(1, len_t)];

% for i=1:gdl
%     figure(i)
%     plot(t, xg(i, :))
% 
% end

x0 = [0; 0; 0; 0];
xd0 = [0; 0; 0; 0];

ZERO = zeros(gdl, gdl);
I = eye(gdl, gdl);
% 
% [t, x_ode] = ode45(@(t, x_ode) ...
%     vibforz_geo(t, x_ode, M, C, K, ZERO, I, L, wl, phase, A), t, [x0; xd0]);

[t, x_ode] = ode45(@(t, x2_ode) ...
    vibforz_geo(t, x2_ode, M, C, K, ZERO, I, L, wl, phase, A), t, [x0; xd0]);

x_ode = x_ode';

for i=1:gdl
    figure(i)
    plot(t, x_ode(i, :))
    hold on

end

clear x_ode
clear xg

close all

%% Método de diferencia central para respuesta ante carga general

clc

amp = 0.01;
b = 0.2;

% Velocidades
% v=[5 20 60 100];
v = 5;
v = v/3.6;
p = b./(dt*v);

da = (l1 + l2)./(dt*v);

for j=1:length(v)
    
    %% Carga
    P0 = zeros(gdl, len_t);
    
    P0(1, 2:round(p(j) + 1)) = triang(round(p(j)))*amp;
    P0(2, (round(1 + da(j)):round((da(j) + p(j))))) = triang(round(p(j)))*amp;

    figure(j)
    plot(t, P0(1, :))
    hold on 
    plot(t, P0(2, :))
    
    P0 = L*P0;  % En [N], lo anterior era solo amplitud
    
    %% Respuesta con diferencia central

    x0 = [0; 0; 0; 0];
    xd0 = [0; 0; 0; 0];
    
    x_difcentral = zeros(4, len_t);
    x_difcentral(:, 1) = x0;

    xd_difcentral = xd0;

    for i=2:len_t
        x_difcentral(:, i) = 0.5*dt^2*(M\P0(:, i) -M\C*xd_difcentral ...
            -M\K*x_difcentral(:, i-1)) + x_difcentral(:, i-1) ... 
            + dt*xd_difcentral;

        xd_difcentral = (2/dt)*(x_difcentral(:, i) - x_difcentral(:, i-1)) ...
            - xd_difcentral;
    end

    for i=1:gdl
        figure(1 + i)
        plot(t, x_difcentral(i, :))
    end
    
end

% v5=zeros(2,len_t);
% v5(1,1:p(1))=triang(p(1)).*0.1;
% v5(2,da(1):(da(1)+p(1)-1))=triang(p(1)).*0.1;
% 
% v20=zeros(2,len_t);
% v20(1,1:p(2))=triang(p(2)).*0.1;
% v20(2,da(2):(da(2)+p(2)-1))=triang(p(2)).*0.1;
% 
% v60=zeros(2,len_t);
% v60(1,1:p(3))=triang(p(3)).*0.1;
% v60(2,da(3):(da(3)+p(3)-1))=triang(p(3)).*0.1;
% 
% v100=zeros(2,len_t);
% v100(1,1:p(4))=triang(p(4)).*0.1;
% v100(2,da(4):(da(4)+p(4)-1))=triang(p(4)).*0.1;
% 



% for i=1:length(v)
%         
%     P0(1, 1:p(i), i) = triang(p(i))*amp;
%     P0(2, (da(i):(da(i)+p(i)-1)), i) = triang(p(i))*amp;
%     
% end
% 
% for i=1:length(v)
%     
%     figure(8+i)
%     plot(t, P0(1, i, :))
%     
% end
    

% figure (1) 
% hold on
% grid on
% ax1=subplot(2,2,1)
% plot(ax1,t,v5,'r');
% 
% ax2=subplot(2,2,2)
% plot(ax2,t,v20,'b');
% ax3=subplot(2,2,3)
% plot(ax3,t,v60,'y');
% ax4=subplot(2,2,4)
% plot(ax4,t,v100,'g');


% x_difcentral = zeros(4, len_t);
% x_difcentral(:, 1) = x0;
% 
% xd_difcentral = xd0;
% 
% for i=2:len_t
%     
%     x_difcentral(:, i) = 0.5*dt^2*(M\P0(:, i) -M\C*xd_difcentral ...
%         -M\K*x_difcentral(:, i-1)) + x_difcentral(:, i-1) ... 
%         + dt*xd_difcentral;
%     
%     xd_difcentral = (2/dt)*(x_difcentral(:, i) - x_difcentral(:, i-1)) ...
%         - xd_difcentral;
%     
% end
% 
% for i=1:gdl
%     figure(i)
%     plot(t, x_difcentral(i, :))
% 
% end
% 
% for i=1:gdl
%     figure(4+i)
%     plot(t, x_difcentral(i, :)-x2_ode(i, :))
% end




%% Función para el ode45 en coordenadas geométricas
%

function x_forz = vibforz_geo(t, x, M, C, K, ZERO, I, L, wl, phase, amp)

    A = [ZERO, I;
        -M\K, -M\C];
    
    B = [zeros(4, 1);
        M\(L*[amp*sin(wl*t);
                amp*sin(wl*t - phase);
                0;
                0])];
    
    x_forz = A*x + B;

end
