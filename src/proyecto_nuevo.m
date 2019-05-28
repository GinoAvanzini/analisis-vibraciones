
clc
clear all
close all

%%
% Datos del sistema

m1 = 25;    % Kg
m2 = m1;    % Kg
m3 = 2000;  % Kg
J = 500;    % Kg*m^2

k1 = 80;    % KN/m
k2 = 50;    % KN/m
k3 = 10;    % KN/m
k4 = k3;    % KN/m

c1 = 4000;  % N*s/m
c2 = 4000;  % N*s/m

l1 = 2;     % m
l2 = 1;     % m

gdl = 4;

dt = 0.005;
t = 0:dt:100;

% Condiciones iniciales en coordenadas geométricas
x0 = [0.002;
    -0.001;
    0;
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

L = [k2;   % Matriz de participación de carga. Pn = L*P(t)
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

%for i = 1:gdl
%     
%     figure(i)
%     plot(t, x(i, :))
%     
% end

%% 
% Respuesta a carga periódica armónica

A = 0.04;  
lambda = 0.2;

v = 40*1000/3600;       % 40 km/h en m/s
wl = 2*pi*v/lambda;     % frec de la carga, función de la vel del auto

beta_n = wl./(Wn_vec.^0.5);
yf = zeros(4, length(t));


P0f = FI'*(A*L);        % Carga generalizada
                              %Hacer una carga con K*x*sen(omega*t)
D = ((1-beta_n.^2).^2 + (2.*zeta_n.*beta_n).^2).^-0.5;
theta = atan(2.*zeta_n.*beta_n./(1-beta_n.^2));


for i=1:gdl
    
   yf(i, :) = (P0f(i)*D(i)/Kn(i)*sin(wl*t-theta(i)));
    
end
P = zeros(4, length(t));
P(:,1)=P0f;
x(:,1)=x0;    %primer columna
%Para Mn(1)
for i=2:length(t)
  x(:,i)=(dt^2/(2*Mn(1)))*((A*L*x(i-1).*sin(wl.*t(i)))-C*xd0-K*x(:,i-1))+x(:,i-1)+dt.*xd0; %Hay que variar el P0f
  xd0=(2/dt)*(x(:,i)-x(:,i-1))-xd0;  
endfor
%Mn(2)
##for i=2:length(t)
##  x(:,i)=(dt^2/(2*Mn(2)))*((A*L*x(i-1).*sin(wl.*t(i)))-C*xd0-K*x(:,i-1))+x(:,i-1)+dt.*xd0; %Hay que variar el P0f
##  xd0=(2/dt)*(x(:,i)-x(:,i-1))-xd0;  
##endfor
%Mn(3)
##for i=2:length(t)
##  x(:,i)=(dt^2/(2*Mn(3)))*((A*L*x(i-1).*sin(wl.*t(i)))-C*xd0-K*x(:,i-1))+x(:,i-1)+dt.*xd0; %Hay que variar el P0f
##  xd0=(2/dt)*(x(:,i)-x(:,i-1))-xd0;  
##endfor
%Mn(4)
##for i=2:length(t)
##  x(:,i)=(dt^2/(2*Mn(4)))*((A*L*x(i-1).*sin(wl.*t(i)))-C*xd0-K*x(:,i-1))+x(:,i-1)+dt.*xd0; %Hay que variar el P0f
##  xd0=(2/dt)*(x(:,i)-x(:,i-1))-xd0;  
##endfor
%Cambiar porque no está bien aplicado 
for  i=1:gdl  
     figure(i)
    plot(t, x(i, :))
    
end
