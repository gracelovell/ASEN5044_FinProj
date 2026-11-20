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

% Given initial state in Earth-centered coordinates
X_0 = r0; % [km] 
Xdot_0 = 0; % [km/s] 
Y_0 = 0; % [km]
Ydot_0 = r0*sqrt(mu/r0^3); % [km/s]
% Pack into initial state vector
x0 = [X_0 Xdot_0 Y_0 Ydot_0]';
n = length(x0); % number of states

% Define time for simulation
delT = 10; % [s] time step
T = 2*pi*sqrt(r0^3/mu); % [s] orbit period
t = 10:delT:14000; % time vector
N = length(t);

% Initialize perturbationn history matrix
dxHist = zeros(n,N);
%%%%%% Choose reasonable initial perturbation?
dx0 = [.1 .1 .005 .005]';
% Initialize as first perturbation in loop
dx_k = dx0;

for k = 1:N
    % Compute linearized DT dynamics matrix for this time step
    alpha = w0*t(k);
    F_k = [1 delT 0 0;
          delT*mu/r0^3*(3*cos(alpha)^2-1)   1   3*delT*mu/r0^3*cos(alpha)*sin(alpha)  0;
          0 0 1 delT;
          3*delT*mu/r0^3*cos(alpha)*sin(alpha)  0  delT*mu/r0^3*(3*sin(alpha)^2-1)  1];
    
    % Compute perturbation at time k using previous perturbation
    dx_k = F_k*dx_k; % No control input perturbations
    
    % Log perturbation in history matrix
    dxHist(:,k) = dx_k;
end

% Generate nominal trajectory
xnomHist = zeros(n,N);
for k = 1:N
    alpha = w0*t(k);
    xnomHist(:,k) = [r0*cos(alpha); -r0*w0*sin(alpha); ...
                     r0*sin(alpha); r0*w0*cos(alpha)];
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

% Model: dy_k = H_k*dx_k

% H_k = a stack C matrices for visible stations at that time tk
% Every time tk has p(tk) = 3*Ns(tk) measurements where Ns is the number of
% visible stations at that time

% A station is visible and generates a 3x1 measurement vector if
% phi_i is within [-pi/2+theta_i, pi/2+theta_i]
% where theta_i(t) = atan2(Y_i(t),X_i(t))

% 12 stations within constant positions on Earth, but rotating because
% Earth is rotating

% At every time t, need to compute new positions for each station and
% evaluate if they are visible, and if so generate a measurement (need to
% build the C matrix for that station at that time and stack with the rest
% to multiply by the state at that time (and its all actually
% perturbations that we add to the nominal, so must generate a nominal 
% solution too))

% Compute initial position for each station
theta0 = zeros(12,1);
Xs_0 = zeros(12,1); % station X position
Ys_0 = zeros(12,1); % station Y position
for i = 1:12
    theta0(i) = (i-1)*pi/6;
    Xs_0(i) = R_E*cos(theta0(i));
    Ys_0(i) = R_E*sin(theta0(i));
end

% Initialize measurement matrices
% A maximum of 2 stations are visible at a time, so need 2*3 rows
dyHist = NaN(2*3,N);
ynomHist = NaN(2*3,N);
yHist = NaN(2*3,N);

% Outer loop: loop through all instants of time
for k = 1:N
    
    % Extract states for this instant
    X_k = xLinHist(1,k);
    Xdot_k = xLinHist(2,k);
    Y_k = xLinHist(3,k);
    Ydot_k = xLinHist(4,k);

    % Extract nominal states for this instant
    Xnom_k = xnomHist(1,k);
    XdotNom_k = xnomHist(2,k);
    Ynom_k = xnomHist(3,k);
    YdotNom_k = xnomHist(4,k);
    H_k = [];
    ynom_k = [];

    % Inner loop: loop through each station
    for i = 1:12

        % Compute position and velocity for each station at this instant
        Xs_i = R_E*cos(wE*t(k) + theta0(i));
        Xsdot_i = -R_E*sin(wE*t(k) + theta0(i))*wE;
        Ys_i = R_E*sin(wE*t(k) + theta0(i));
        Ysdot_i = R_E*cos(wE*t(k) + theta0(i))*wE;
        theta_i = atan2(Ys_i,Xs_i);

        % Compute phi to check station visibility
        phi_i = atan2(Ynom_k - Ys_i,Xnom_k - Xs_i);

        % Wrap the relative angle (phi_i - theta_i) into [-pi, pi]
        relAngle = atan2(sin(phi_i - theta_i), cos(phi_i - theta_i));

        % Check visibility
        if relAngle >= -pi/2 && relAngle <= pi/2
            % If station is visible, compute rhoNom and rhodotNom to build
            % C matrix
            rhoNom_i = sqrt((Xnom_k - Xs_i)^2 + (Ynom_k - Ys_i)^2);
            rhodotNom_i = ((Xnom_k - Xs_i)*(XdotNom_k - Xsdot_i) ...
                + (Ynom_k - Ys_i)*(YdotNom_k - Ysdot_i))/rhoNom_i;

            % Store nominal measurement
            ynom_i = [rhoNom_i rhodotNom_i phi_i]';
            ynom_k = [ynom_k; ynom_i];
            
            % Build C matrix
            C11 = (Xnom_k - Xs_i)/rhoNom_i;
            C13 = (Ynom_k - Ys_i)/rhoNom_i;
            C21 = (XdotNom_k - Xsdot_i)/rhoNom_i ...
                   - (rhodotNom_i*(Xnom_k - Xs_i))/rhoNom_i^2;
            C22 = C11;
            C23 = (YdotNom_k - Ysdot_i)/rhoNom_i ...
                   - (rhodotNom_i*(Ynom_k - Ys_i))/rhoNom_i^2;
            C24 = C13;
            C31 = -(Ynom_k - Ys_i)/rhoNom_i^2;
            C33 = (Xnom_k - Xs_i)/rhoNom_i^2;
            C_i = [C11 0 C13 0; C21 C22 C23 C24; C31 0 C33 0];
            
            % Stack C matrices into H_k matrix
            H_k = [H_k; C_i];
            % H_k will have 3*Ns where Ns is the number of visible stations
        end
            % if  station is not visible, no C matrix, no measurement
    end

    % Compute measurement perturbation for this instant
    dy_k = H_k*dxHist(:,k); % dy_k will be 3*Ns x 1
    % Store measurement perturbation 
    dyHist(1:length(dy_k),k) = dy_k;

    % Compute total measurement for this instant
    y_k = dy_k + ynom_k;
    % Keep phi measurements between -pi and pi
    y_k(3) = atan2(sin(y_k(3)), cos(y_k(3)));
    if length(y_k) == 6
        y_k(6) = atan2(sin(y_k(6)), cos(y_k(6)));
    end
    % Store nominal and total measurements
    ynomHist(1:length(ynom_k),k) = ynom_k;
    yHist(1:length(y_k),k) = y_k;
end

% Plot simulated measurements
figure
% Plot rho measurements with nominal measurements
subplot(3,1,1)
plot(t,yHist(1,:),'x','Color','blue')
hold on
plot(t,yHist(4,:),'x','Color','blue')
xlabel('Time (s)')
ylabel('\rho Measurements [UNITS]')
title('Simulated Linearized \rho Measurements')
% Plot rhodot measurements with nominal measurements
subplot(3,1,2)
plot(t,yHist(2,:),'o','Color','blue')
hold on
plot(t,yHist(5,:),'o','Color','blue')
xlabel('Time (s)')
ylabel('\rho dot Measurements [UNITS]')
title('Simulated Linearized \rho dot Measurements')
% Plot phi measurements with nominal measurements
subplot(3,1,3)
plot(t,yHist(3,:),'s','Color','blue')
hold on
plot(t,yHist(6,:),'s','Color','blue')
xlabel('Time (s)')
ylabel('\phi Measurements [UNITS]')
title('Simulated Linearized \phi Measurements')







