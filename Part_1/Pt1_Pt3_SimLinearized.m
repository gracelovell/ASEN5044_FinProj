% Pt 1: take given nonlinear dynamics and create CT linearized model 
% matrices. This means saying the state is approximately equal
% to the nominal state plus perturbations, and state derivative is
% equal to the nonlinear dynamics evaluated at the nominal state +
% perturbations + the jacobian of the nonlinear dynamics evaluated at the
% nominal state times the perturbations, plus higher order terms that we
% ignore. then we find that the derivative of perturbations is equalt to
% those jacobians times the perturbations. So we can compute the
% change in perturbations instead of the change in state. Then all of the
% jacobians are computed to find A, B, C, D, the continuous time linearized
% model matrices

% Pt 2: Finding the nominal state and evaluating the CT matrices here. The
% initial state is given, and we use 0 input since it is a natural circular
% orbit. The nominal trajectory needs to be computed in X and Y
% coordinates, so it is not a constant state (r and v would be constant in
% rotating rth frame). So the nominal state is a function of time, and we
% can evaluate the matrices here.
% Then, we convert to DT. This means we got from dxdot = Adx + Bdu to dx(k+1) =
% Fdx(k) + Gdu(k). The new matrices are approximated with Euler integration

% Pt 3: Code up the DT model matrices and simulate the linerized dynamics.
% Also simulate the nonlinear dynamics.

%% Simulate Linearized DT Dynamics

% Given constant parameters
mu = 398600; % [km^3/s^2] gravitational parameter for Earth
r0 = 6678; % [km] nominal orbit radius
R_E = 6378; % [km] Earth's radius
wE = 2*pi/86400; % [rad/s] Earth's angular rotation rate
w0 = sqrt(mu/r0^3); % [rad/s] spacecraft's nominal angular rate

% Given initial state
x1_0 = r0; % [km] X in Earth-centered coordinates
x2_0 = 0; % [km/s] Xdot in Earth-centered coordinates
x3_0 = 0; % [km] Y in Earth-centered coordinates
x4_0 = r0*sqrt(mu/r0^3); % [km/s] Ydot in Earth-centered coordinates
% Pack into initial state vector
x0 = [x1_0 x2_0 x3_0 x4_0]';
n = length(x0); % number of states

% Define time for simulation
delT = 10; % [s] time step
T = 2*pi*sqrt(r0^3/mu); % [s] orbit period
N = floor(T/delT) + 1; % number of simulation steps
t = (0:N-1)*delT; % time vector

% Initialize perturbationn history matrix
dxHist = zeros(n,N);
%%%%%% Choose reasonable initial perturbation?
dx0 = [.1 .1 .01 .01]';
dxHist(:,1) = dx0;
% Initialize as first perturbation in loop
dx_k = dx0;

for k = 2:N
    % Compute linearized DT dynamics matrix for this time step
    theta = w0*t(k);
    F_k = [1 delT 0 0;
          delT*mu/r0^3*(3*cos(theta)^2-1)   1   3*delT*mu/r0^3*cos(theta)*sin(theta)  0;
          0 0 1 delT;
          3*delT*mu/r0^3*cos(theta)*sin(theta)  0  delT*mu/r0^3*(3*sin(theta)^2-1)  1];
    
    % Compute perturbation at time k using previous perturbation
    dx_k = F_k*dx_k; % No control input perturbations
    
    % Log perturbation in history matrix
    dxHist(:,k) = dx_k;
end

% Generate nominal trajectory
xnomHist = zeros(n,N);
for k = 1:N
    theta = w0*t(k);
    xnomHist(:,k) = [r0*cos(theta); -r0*w0*sin(theta); ...
                     r0*sin(theta); r0*w0*cos(theta)];
end

% Add perturbations to nominal trajectory to generate linearized state
xLinHist = xnomHist + dxHist;

% Plot linearized state history
figure
subplot(2,2,1)
plot(t,xLinHist(1,:))
hold on
plot(t,xnomHist(1,:))
xlabel('Time (s)')
ylabel('X (km)')
title('Evolution Linearized X Position')
legend('Linearized State','Nominal State')

subplot(2,2,2)
plot(t,xLinHist(2,:))
hold on
plot(t,xnomHist(2,:))
xlabel('Time (s)')
ylabel('Xdot (km/s)')
title('Evolution of Linearized X Velocity')
legend('Linearized State','Nominal State')

subplot(2,2,3)
plot(t,xLinHist(3,:))
hold on
plot(t,xnomHist(3,:))
xlabel('Time (s)')
ylabel('Y (km)')
title('Evolution of Linearized Y Position')
legend('Linearized State','Nominal State')

subplot(2,2,4)
plot(t,xLinHist(4,:))
hold on
plot(t,xnomHist(4,:))
xlabel('Time (s)')
ylabel('Ydot (km/s)')
title('Evolution of Linearized Y Velocity')
legend('Linearized State','Nominal State')

%% Simulate Linearized DT Measurements