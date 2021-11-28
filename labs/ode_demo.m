%% Initial Value Problems using ODE45

%% Scalar equation
% dy/dt = sin( (t+y)^2 ) (nonlinear)

f = @(t, y) sin( (t+y).^2 );
tspan = [0, 4];
% tspan = linspace(0, 4, 101);
y0 = 0;
opts = odeset('AbsTol', 1e-8, 'RelTol', 1e-5);
[t, y] = ode45(f, tspan, y0, opts);

plot(t, y), grid on
xlabel('t')
ylabel('y')

%% Direction Fields
tv = linspace(0, 4, 21);
yv = linspace(-2, 1, 21);
[T, Y] = meshgrid(tv, yv);
DY = f(T, Y);
DT = ones(size(DY));
quiver(T, Y, DT, DY);
axis equal, axis tight
set(gca, 'TickDir', 'Out')
title('Direction Field')
xlabel('t')
ylabel('y')

hold on
plot(t, y, 'r')


%% System of ODEs

%% Spring-Mass-Damper system
clear, clf
omega0 = 1;                             % natural freq
zeta = 0.2;                             % damping
A = [0 1; -omega0^2 -2*zeta*omega0];
f = @(t, y) A*y;
tspan = [0 30];
y0 = [3; 0];                            % y1(0) = x = 3, y2(0) = x' = 0

% ODE45
[t, y] = ode45(f, tspan, y0);
plot(t, y)

%%
% Euler, Euler-Trapezoidal, RK4
nstep = 500;
t = linspace(tspan(1), tspan(2), nstep+1);
h = t(2)-t(1);
y1 = zeros(2, nstep+1);                 % 1st row: y1, 2nd row: y2
y1(:,1) = y0;
y2 = y1;

for j = 1:nstep
    y1(:,j+1) = y1(:,j) + f(t(j), y1(:,j))*h; % Forward Euler
    
    s1 = f(t(j), y2(:,j));
    s2 = f(t(j)+h, y2(:,j) + s1*h);
    y2(:,j+1) = y2(:,j) + (s1 + s2)*h/2; % Euler-Trap.
end
[~, y4] = ode45(f, t, y0);

figure(1)
plot(t, y1(1,:), t, y2(1,:), t, y4(:,1))
legend('1st-order', '2nd-order', '4th-order');
title('x(t) using three methods')

% figure(2)
% plot(y1(1,:), y1(2,:));
% xlabel('x'), ylabel('v')
% title('phase diagram')


