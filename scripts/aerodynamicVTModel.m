% --- Modelo aerodinámico del estabilizador vertical con entrada de control ---

% Datos del perfil y geometría del estabilizador vertical (VT)
perfilVT = 'NACA0010'; % Perfil del estabilizador vertical
cuerdaRaiz = 2.8; % Cuerda en la raíz en metros
cuerdaPunta = 1.0; % Cuerda en la punta en metros
alturaVT = 5.5; % Altura total del estabilizador vertical en metros
span= alturaVT*2;
segments = 20; % Número de segmentos en el estabilizador vertical
camberPorcentage = 0; % Sin camber (plano) para el estabilizador vertical

% Parámetros geométricos del VT
AR_VT = (span^2) / ((cuerdaRaiz + cuerdaPunta) / 2 * alturaVT); % Relación de aspecto
LocalChord_VT = linspace(cuerdaRaiz, cuerdaPunta, segments); % Cuerda local de cada segmento
SurfaceLocal_VT = zeros(1, segments); % Inicialización de superficies locales de cada segmento
e = 1.78 * (1 - 0.045 * AR_VT^0.68) - 0.64; % Factor de eficiencia de Oswald estimado para el estabilizador vertical

% Cálculo de superficie total y local usando método de trapecio
h = alturaVT / segments; % Altura de cada segmento
for i = 1:segments
    SurfaceLocal_VT(i) = ((LocalChord_VT(i) + LocalChord_VT(min(i+1, segments))) / 2) * h;
end
SurfaceTotal_VT = sum(SurfaceLocal_VT); % Superficie total del estabilizador vertical

% Velocidades de interés
velocidades = [10, 80, 210, 270]; % Velocidades en m/s
AoA_general = atan(10/100); % Ángulo de ataque general en radianes

% Condiciones atmosféricas usando la función atmosmodel
altitud = 10000; % Altitud de crucero en metros
temperatura = 15; % Temperatura en °C a altitud de crucero
presion_local = 0; % Presión local en metros para referencia de nivel del mar
[T, P, rho, mu] = atmosmodel(altitud, temperatura, presion_local);

% Parámetros aerodinámicos
CorrectedSlope = (2 * pi) / (1 + (2 * pi) / (pi * AR_VT)); % Pendiente corregida del estabilizador

% Entrada de control: Deflexión del estabilizador vertical
deflexion_estabilizador = 5; % Deflexión del estabilizador en grados
deflexion_estabilizador_rad = deflexion_estabilizador * pi / 180; % Conversión a radianes

% --- Inicialización de vectores de resultados ---
Reynolds_VT = zeros(1, segments);
Cl_local_VT = zeros(1, segments);
Cd_friction_VT = zeros(1, segments);
D_friction_VT = zeros(1, segments);
Cd_induced_VT = zeros(1, segments);
D_induced_VT = zeros(1, segments);
L_local_VT = zeros(1, segments);

% --- Calcular coeficientes de sustentación, fricción y arrastre inducido para cada segmento ---
for i = 1:segments
    % Calcular número de Reynolds para cada segmento
    Reynolds_VT(i) = rho * velocidades(2) * LocalChord_VT(i) / mu;
    
    % Coeficiente de fricción
    if Reynolds_VT(i) < 3e6
        Cd_friction_VT(i) = 1.328 / sqrt(Reynolds_VT(i)); % Flujo laminar
    else
        Cd_friction_VT(i) = 0.074 / Reynolds_VT(i)^(1/5); % Flujo turbulento
    end
    
    % Fuerza de fricción
    D_friction_VT(i) = 0.5 * rho * velocidades(2)^2 * SurfaceLocal_VT(i) * Cd_friction_VT(i);
    
    % Ángulo de ataque local ajustado por la deflexión del estabilizador
    AoALocal = AoA_general + deflexion_estabilizador_rad; % Ajuste del ángulo de ataque debido a la deflexión
    
    % Coeficiente de sustentación local
    Cl_local_VT(i) = CorrectedSlope * AoALocal;
    
    % Fuerza de sustentación local
    L_local_VT(i) = 0.5 * rho * velocidades(2)^2 * SurfaceLocal_VT(i) * Cl_local_VT(i);
    
    % Coeficiente de arrastre inducido
    Cd_induced_VT(i) = Cl_local_VT(i)^2 / (pi * AR_VT * e);
    
    % Fuerza de arrastre inducido
    D_induced_VT(i) = 0.5 * rho * velocidades(2)^2 * SurfaceLocal_VT(i) * Cd_induced_VT(i);
end

% --- Sumar sustentación total y arrastre total en el estabilizador vertical ---
LiftTotal_VT = sum(L_local_VT);
DragTotal_friction_VT = sum(D_friction_VT);
DragTotal_induced_VT = sum(D_induced_VT);
DragTotal_VT = DragTotal_friction_VT + DragTotal_induced_VT;

% --- Generación de la Distribución de Sustentación y Arrastre para el Estabilizador Vertical Completo ---
% No es necesario duplicar como en el HT, pues el VT es simétrico sobre el eje vertical
y_VT = linspace(0, alturaVT, segments); % Posición a lo largo del estabilizador vertical

% --- Gráfica de sustentación en el estabilizador vertical completo ---
figure;
plot(y_VT, L_local_VT, 'b-', 'DisplayName', 'Sustentación Estabilizador Vertical');
xlabel('Posición a lo Largo del Estabilizador Vertical (m)');
ylabel('Sustentación (N)');
title(['Distribución de Sustentación en el Estabilizador Vertical con Deflexión de ', num2str(deflexion_estabilizador), '°']);
grid on;
legend;

% --- Gráfica de arrastre total en el estabilizador vertical completo ---
figure;
plot(y_VT, D_friction_VT + D_induced_VT, 'r-', 'DisplayName', 'Arrastre Total Estabilizador Vertical');
xlabel('Posición a lo Largo del Estabilizador Vertical (m)');
ylabel('Arrastre Total (N)');
title(['Distribución de Arrastre Total en el Estabilizador Vertical con Deflexión de ', num2str(deflexion_estabilizador), '°']);
grid on;
legend;

% --- Imprimir resultados ---
fprintf('Número de Reynolds en cada segmento: %.2e\n', Reynolds_VT);
fprintf('Sustentación total del estabilizador vertical con deflexión de %.1f°: %.2f N\n', deflexion_estabilizador, LiftTotal_VT);
fprintf('Arrastre total (fricción) del estabilizador vertical: %.2f N\n', DragTotal_friction_VT);
fprintf('Arrastre total (inducido) del estabilizador vertical: %.2f N\n', DragTotal_induced_VT);
fprintf('Arrastre total del estabilizador vertical: %.2f N\n', DragTotal_VT);
