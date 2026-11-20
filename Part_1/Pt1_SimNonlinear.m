clc
clear
close all

% givens
mu = 398600; %km^3/s
r0 = sqrt(6678^2 + 0^2);

x0 = [6678; 0; 0; r0*sqrt(mu/r0^3)];

% simulate nonlinear dynamics
tspan = [0 14000];
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
[t, x_vec] = ode45(@(t, x) nonlin_dynam(mu, x), tspan, x0, options);

% plot results
figure(1)

% X coordinate
subplot(4, 1, 1)
plot(t, x_vec(:, 1), LineWidth=1.5)
xlabel("Time (secs)")
ylabel("X (km)")

% Xdot coordinate
subplot(4, 1, 2)
plot(t, x_vec(:, 2), LineWidth=1.5)
xlabel("Time (secs)")
ylabel("Xdot (km/s)")

subplot(4, 1, 3)
plot(t, x_vec(:, 3), LineWidth=1.5)
xlabel("Time (secs)")
ylabel("Y (km)")

subplot(4, 1, 4)
plot(t, x_vec(:, 4), LineWidth=1.5)
xlabel("Time (secs)")
ylabel("Ydot (km/s)")


% function for ode45 to integrate
function[x_dot] = nonlin_dynam(mu, x)

x_dot = [x(2)
    -mu*x(1)/((sqrt(x(1)^2 + x(3)^2))^3)
    x(4)
    -mu*x(3)/((sqrt(x(1)^2 + x(3)^2))^3)];

end
