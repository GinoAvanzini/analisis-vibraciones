%%
Respuesta en vibraciones libres analíticamente

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
Respuesta a carga periódica armónica de forma analítica


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
Respuesta a carga armónica con ode45


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
Respuesta a carga impulsiva mediante la Integral de Duhamel

a = 0.05;
b = 0.04;

v = 20*1000/3600;

% pico = ((a/2)/v)

% tri = triangularPulse(0, (a/2)/v, 0.8, t);

% figure(8)
% plot(t, tri)

P = k2*A.*sin(wl*t);


y_duha = zeros(4, len_t);

for i=1:gdl
   
   yc = P.*cos(Wn_vec(i)^0.5*t);
   ys = P.*sin(Wn_vec(i)^0.5*t);

   
   A_d = zeros(1, len_t);
   B_d = zeros(1, len_t);

   decr = exp(-zeta_n(i)*Wn_vec(i)^0.5*dt);
   
   for j = 2:len_t
       
       A_d(j) = A_d(j-1)*decr + ...
           (dt/(2*Mn(i)*Wn_vec(i)^0.5))*(yc(j-1)*decr + yc(j));
       B_d(j) = B_d(j-1)*decr + ...
           (dt/(2*Mn(i)*Wn_vec(i)^0.5))*(ys(j-1)*decr + ys(j));
       
   end
   
   aa = A_d.*sin(Wnd(i)*t);
   bb = B_d.*cos(Wnd(i)*t);
   y_duha(i, :) = aa - bb;

end

for i=1:gdl
   figure(i)
   plot(t, y_duha(i, :))
end




%% Función para el ode45 en coordenadas modales


function yd = vibforz(t, y, z, w, P0m, wL)

    A = [0, 1;
        -w^2, -2*z*w];
    L = [0;
        P0m*sin(wL*t)];
    yd = A*y + L;

end



