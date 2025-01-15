% --- Modelo aerodinámico del estabilizador horizontal con entrada de control ---
function [LiftTotalHT] = aerodynamicHTModel(deflexion_estabilizador)

% Datos del perfil y geometría del estabilizador horizontal (HT)
perfilHT = 'NACA0012'; % Perfil del estabilizador horizontal
cuerdaRaiz = 3.1; % Cuerda en la raíz en metros
cuerdaPunta = 1.3; % Cuerda en la punta en metros
envergaduraHT = 8.0; % Envergadura total del estabilizador horizontal en metros
segments = 20; % Número de segmentos en el estabilizador horizontal
camberPorcentage = 0; % Sin camber (plano) para el estabilizador horizontal

% Parámetros geométricos del HT
AR_HT = (envergaduraHT^2) / ((cuerdaRaiz + cuerdaPunta) / 2 * envergaduraHT); % Relación de aspecto
LocalChord_HT = linspace(cuerdaRaiz, cuerdaPunta, segments); % Cuerda local de cada segmento
SurfaceLocal_HT = zeros(1, segments); % Inicialización de superficies locales de cada segmento
e = 1.78 * (1 - 0.045 * AR_HT^0.68) - 0.64; % Factor de eficiencia de Oswald estimado para el estabilizador horizontal

% Cálculo de superficie total y local usando método de trapecio
h = envergaduraHT / (2 * segments); % Ancho de cada segmento
for i = 1:segments
    SurfaceLocal_HT(i) = ((LocalChord_HT(i) + LocalChord_HT(min(i+1, segments))) / 2) * h;
end
SurfaceTotal_HT = sum(SurfaceLocal_HT); % Superficie total del estabilizador horizontal

% Velocidades de interés
velocidades = [10, 80, 210, 270]; % Velocidades en m/s
AoA_general = atan(10/100); % Ángulo de ataque general en radianes

% Condiciones atmosféricas usando la función atmosmodel
altitud = 10000; % Altitud de crucero en metros
temperatura = 15; % Temperatura en °C a altitud de crucero
presion_local = 0; % Presión local en metros para referencia de nivel del mar
[T, P, rho, mu] = atmosmodel(altitud, temperatura, presion_local);

% Parámetros aerodinámicos
CorrectedSlope = (2 * pi) / (1 + (2 * pi) / (pi * AR_HT)); % Pendiente corregida del estabilizador

% Entrada de control: Deflexión del estabilizador horizontal
% Deflexión del estabilizador en grados
deflexion_estabilizador_rad = deflexion_estabilizador * pi / 180; % Conversión a radianes

% --- Inicialización de vectores de resultados ---
Reynolds_HT = zeros(1, segments);
Cl_local_HT = zeros(1, segments);
Cd_friction_HT = zeros(1, segments);
D_friction_HT = zeros(1, segments);
Cd_induced_HT = zeros(1, segments);
D_induced_HT = zeros(1, segments);
L_local_HT = zeros(1, segments);

% --- Calcular coeficientes de sustentación, fricción y arrastre inducido para cada segmento ---
for i = 1:segments
    % Calcular número de Reynolds para cada segmento
    Reynolds_HT(i) = rho * velocidades(2) * LocalChord_HT(i) / mu;
    
    % Coeficiente de fricción
    if Reynolds_HT(i) < 3e6
        Cd_friction_HT(i) = 1.328 / sqrt(Reynolds_HT(i)); % Flujo laminar
    else
        Cd_friction_HT(i) = 0.074 / Reynolds_HT(i)^(1/5); % Flujo turbulento
    end
    
    % Fuerza de fricción
    D_friction_HT(i) = 0.5 * rho * velocidades(2)^2 * SurfaceLocal_HT(i) * Cd_friction_HT(i);
    
    % Ángulo de ataque local ajustado por la deflexión del estabilizador
    AoALocal = AoA_general + deflexion_estabilizador_rad; % Ajuste del ángulo de ataque debido a la deflexión
    
    % Coeficiente de sustentación local
    Cl_local_HT(i) = CorrectedSlope * AoALocal;
    
    % Fuerza de sustentación local
    L_local_HT(i) = 0.5 * rho * velocidades(2)^2 * SurfaceLocal_HT(i) * Cl_local_HT(i);
    
    % Coeficiente de arrastre inducido
    Cd_induced_HT(i) = Cl_local_HT(i)^2 / (pi * AR_HT * e);
    
    % Fuerza de arrastre inducido
    D_induced_HT(i) = 0.5 * rho * velocidades(2)^2 * SurfaceLocal_HT(i) * Cd_induced_HT(i);
end

% --- Sumar sustentación total y arrastre total en el estabilizador horizontal ---
LiftTotal_HT = sum(L_local_HT);
DragTotal_friction_HT = sum(D_friction_HT);
DragTotal_induced_HT = sum(D_induced_HT);
DragTotal_HT = DragTotal_friction_HT + DragTotal_induced_HT;

% --- Generación de la Distribución de Sustentación y Arrastre para el Estabilizador Completo ---
% Duplicamos la distribución de sustentación y arrastre y reflejamos el lado izquierdo
x_right = linspace(0, envergaduraHT / 2, segments); % Posición a lo largo del estabilizador derecho
x_left = -fliplr(x_right); % Posición a lo largo del estabilizador izquierdo
L_local_left = fliplr(L_local_HT); % Sustentación en el lado izquierdo
D_total_left = fliplr(D_friction_HT + D_induced_HT); % Arrastre total en el lado izquierdo

% Combinar los datos para el gráfico completo
x_full = [x_left, x_right]; % Coordenadas x para ambos lados
L_full = [L_local_left, L_local_HT]; % Sustentación en ambos lados
D_full = [D_total_left, D_friction_HT + D_induced_HT]; % Arrastre total en ambos lados

% --- Gráfica de sustentación en el estabilizador horizontal completo ---
figure;
plot(x_full, L_full, 'b-', 'DisplayName', 'Sustentación Estabilizador Horizontal');
xlabel('Posición a lo Largo del Estabilizador (m)');
ylabel('Sustentación (N)');
title(['Distribución de Sustentación en el Estabilizador Horizontal Completo con Deflexión de ', num2str(deflexion_estabilizador), '°']);
grid on;
legend;

% --- Gráfica de arrastre total en el estabilizador horizontal completo ---
figure;
plot(x_full, D_full, 'r-', 'DisplayName', 'Arrastre Total Estabilizador Horizontal');
xlabel('Posición a lo Largo del Estabilizador (m)');
ylabel('Arrastre Total (N)');
title(['Distribución de Arrastre Total en el Estabilizador Horizontal Completo con Deflexión de ', num2str(deflexion_estabilizador), '°']);
grid on;
legend;

% --- Imprimir resultados ---
fprintf('Número de Reynolds en cada segmento: %.2e\n', Reynolds_HT);
fprintf('Sustentación total del estabilizador horizontal con deflexión de %.1f°: %.2f N\n', deflexion_estabilizador, LiftTotal_HT);
fprintf('Arrastre total (fricción) del estabilizador horizontal: %.2f N\n', DragTotal_friction_HT);
fprintf('Arrastre total (inducido) del estabilizador horizontal: %.2f N\n', DragTotal_induced_HT);
fprintf('Arrastre total del estabilizador horizontal: %.2f N\n', DragTotal_HT);
end
