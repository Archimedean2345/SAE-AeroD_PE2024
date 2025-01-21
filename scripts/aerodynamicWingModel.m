function [LiftTotalWing, LiftTotal_right, LiftTotal_left,DragTotalWing, DragTotal_right, DragTotal_left] = aerodynamicWingModel(aleronDeflection)
% --- Modelo aerodinamico de ala ---
% max deflexion of aleron 25deg
% --- Declaration of aerodynamic parameters ---
[ARwing, AoALocal, LocalChord, SurfaceLocal, surfaceWing, tipChordWing, rootChordWing, spanWing, Segments, camberWingPorcentage] = geometricParameters();
[T, P, rho, mu] = atmosmodel(1870,15,0);

% Update camber and corrected slope with new variable names
camber = camberWingPorcentage / 100;
ZeroLiftAngle = -2*camber;
Cd_pressure = 0.02;
ewing= 1.78 * (1 - 0.045 * ARwing^0.68) - 0.64; % Oswald parameter for straight wings and medium size AR
CorrectedSlope = (2*pi)/(1+ (2*pi)/(ewing*pi*ARwing)); % Corrected slope

% Auxiliary variables 
Velocity_u = 20; % Forward velocity in m/s
aleronDeflectionAngle = 20; % Aileron deflection angle in degrees
Velocity_w = 1; % Vertical velocity in m/s

% --- Initialization of Vectors ---
Reynolds = zeros(1, Segments);
D_pressure = 0.5 * rho * Velocity_u^2 * SurfaceLocal * Cd_pressure;
Cd_friction = zeros(1, Segments);
D_friction = zeros(1, Segments);
Cd_induced = zeros(1, Segments);
D_induced = zeros(1, Segments);
D_interference = zeros(1, Segments);
Cl_local = zeros(1, Segments);
L_local = zeros(1, Segments);

% --- General Angle of Attack for all segments ---
AoAGeneral = atan(Velocity_w / Velocity_u) * ones(1, Segments);
AoAGeneral= deg2rad(2)* ones(1, Segments);  % Manually indicated angle of attack

% --- Calculation of Lift and Drag Coefficients for each segment ---
for i = 1:Segments
    % Reynolds number for segment i
    Reynolds(i) = rho * Velocity_u * LocalChord(i) / mu;
    
    % Friction coefficient (laminar/turbulent)
    if Reynolds(i) < 3e6
        Cd_friction(i) = 1.328 / sqrt(Reynolds(i)); % Laminar flow
    else
        Cd_friction(i) = 0.074 / Reynolds(i)^(1/5); % Turbulent flow
    end
    D_friction(i) = 0.5 * rho * Velocity_u^2 * SurfaceLocal(i) * Cd_friction(i);

    % Initial induced drag coefficient (no aileron)
    Cd_induced(i) = Cl_local(i)^2 / (pi * ARwing * ewing);
    D_induced(i) = 0.5 * rho * Velocity_u^2 * SurfaceLocal(i) * Cd_induced(i);

    % Interference drag (no aileron)
    D_interference(i) = 0.2*(D_friction(i) + D_pressure(i) + D_induced(i));

    % Lift coefficient (no aileron)
    Cl_local(i) = CorrectedSlope * (AoAGeneral(i) + AoALocal(i) - ZeroLiftAngle);
end

% --- Lift Coefficient Modification Due to Aileron Deflection ---
[camberAleron, AleronWingInitialLocationPorcentage, AleronWingLastLocationPorcentage] = controlSurfacesModelwing(aleronDeflection);
ZeroLiftAngleRight= -2 *camberAleron;
[camberAleron1] = controlSurfacesModelwing(-aleronDeflection);
ZeroLiftAngleLeft= -2 *camberAleron1;
AleronWingInitialSegment = round(Segments *AleronWingInitialLocationPorcentage /100);
AleronWingLastSegment = round(Segments *AleronWingLastLocationPorcentage /100);

% Initialize Cl and Cd vectors for each wing
Cl_local_right = Cl_local;
Cl_local_left = Cl_local;
Cd_induced_right = zeros(1, Segments);
Cd_induced_left = zeros(1, Segments);

% Opposing aileron deflection on each wing
for i = AleronWingInitialSegment:AleronWingLastSegment
    % Right wing: positive deflection
    Cl_local_right(i) = CorrectedSlope * (AoAGeneral(i) + AoALocal(i) - ZeroLiftAngleRight);
    
    % Left wing: negative deflection
    Cl_local_left(i) = CorrectedSlope * (AoAGeneral(i) + AoALocal(i) - ZeroLiftAngleLeft);
end

% --- Modified induced drag calculation for each wing ---
for i = 1:Segments
    Cd_induced_right(i) = Cl_local_right(i)^2 / (pi*ARwing*ewing);
    Cd_induced_left(i) = Cl_local_left(i)^2 / (pi*ARwing*ewing);
end

% Calculate induced drag forces for each segment on each wing
D_induced_right = 0.5 * rho * Velocity_u^2 * SurfaceLocal .* Cd_induced_right;
D_induced_left = 0.5 * rho * Velocity_u^2 * SurfaceLocal .* Cd_induced_left;

% --- Calculation of Lift Forces on Each Wing ---
L_local_right = 0.5 * rho * Velocity_u^2 .* SurfaceLocal .* Cl_local_right;
L_local_left = 0.5 * rho * Velocity_u^2 .* SurfaceLocal .* Cl_local_left;

% --- Fuerzas sustentadoras y de resistencia ---
LiftTotal_right = sum(L_local_right);
LiftTotal_left = sum(L_local_left);
DragTotal_right = sum(D_induced_right) + sum(D_friction) + sum(D_pressure);
DragTotal_left = D_induced_left + sum(D_friction) + sum(D_pressure);
LiftTotalWing = LiftTotal_left + LiftTotal_right;
DragTotalWing = DragTotal_right + DragTotal_left;

% Compute center positions of lift and drag
x_right = linspace(0, spanWing / 2, Segments);  % Position along right wing
x_left = -linspace(0, spanWing / 2, Segments);  % Position along left wing

LiftCenter_right = sum(x_right .* L_local_right) / LiftTotal_right;
LiftCenter_left = sum(x_left .* L_local_left) / LiftTotal_left;
DragCenter_right = sum(x_right .* D_induced_right) / sum(D_induced_right);
DragCenter_left = sum(x_left .* D_induced_left) / sum(D_induced_left);

% --- Sustentacion ---
figure;
hold on;
plot(x_right, L_local_right, 'b-', 'DisplayName', 'Right Wing Lift');
plot(x_left, L_local_left, 'r-', 'DisplayName', 'Left Wing Lift');
xline(LiftCenter_right, 'b--', ['Right Lift Center: ', num2str(LiftCenter_right, '%.2f'), ' m']);
xline(LiftCenter_left, 'r--', ['Left Lift Center: ', num2str(LiftCenter_left, '%.2f'), ' m']);
xlabel('Position Along Wing (m)'); 
ylabel('Lift (N)');
title('Lift Distribution Across the Wing');
grid on;
hold off;

% --- Resistencia al avance ---
figure;
hold on;
plot(x_right, D_induced_right, 'b-', 'DisplayName', 'Right Wing Induced Drag');
plot(x_left, D_induced_left, 'r-', 'DisplayName', 'Left Wing Induced Drag');
xline(DragCenter_right, 'b--', ['Right Drag Center: ', num2str(DragCenter_right, '%.2f'), ' m']);
xline(DragCenter_left, 'r--', ['Left Drag Center: ', num2str(DragCenter_left, '%.2f'), ' m']);
xlabel('Position Along Wing (m)');
ylabel('Drag (N)');
title('Induced Drag Distribution Across the Wing');
grid on;
hold off;
end
