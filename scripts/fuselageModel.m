function fuselageModel(numPuntos)
    % fuselajeTrainerJet - Genera un modelo de fuselaje estilo wireframe para un avión trainer jet.
    % numPuntos: Número de puntos para discretizar cada elipse.

    % Longitudes de las secciones en el eje Z (longitudinal)
    seccionNariz = 2;            % Longitud de la nariz (0 a 2 metros).
    seccionCabina = [2, 3.8];    % Cabina de la tripulación (2 a 3.8 metros).
    seccionTransicion = [3.8, 6]; % Zona de transición (3.8 a 6 metros).
    seccionEmpenaje = [6, 10];   % Empenaje trasero (6 a 10 metros).

    % Parámetros de elipses en cada sección (a: semi-eje mayor, b: semi-eje menor)
    radioNariz = [0.2, 0.1];     % Radio de la elipse de la nariz (punta hacia una elipse)
    radioCabina = [0.7, 0.5];    % Radio de la elipse de la cabina
    radioTransicion = [0.5, 0.35]; % Radio de la elipse de transición
    radioEmpenaje = [0.2, 0.1];  % Radio final en la cola del avión

    % Posición vertical de cada elipse para simular elevación en la cola
    elevacionCabina = 0;
    elevacionTransicion = 0;
    elevacionEmpenaje = 0.3;     % Desplazamiento hacia arriba en la cola

    % Configuración de la figura
    figure;
    hold on;
    axis equal;
    xlabel('Ancho del fuselaje (m)');
    ylabel('Altura del fuselaje (m)');
    zlabel('Longitud del fuselaje (m)');
    title('Modelo Wireframe de un Fuselaje Trainer Jet');
    grid on;

    % 1. Nariz puntiaguda con elipse de base
    [x, y, z] = crearNariz(seccionNariz, radioNariz, numPuntos);
    plot3(x, y, z, 'k'); % Dibuja la nariz en estilo wireframe

    % 2. Cabina
    posicionesCabina = linspace(seccionCabina(1), seccionCabina(2), 5);
    elipsesCabina = [];
    for pos = posicionesCabina
        [x, y, z] = elipse3D(radioCabina(1), radioCabina(2), pos, elevacionCabina, numPuntos);
        elipsesCabina = [elipsesCabina; x', y', z'];
        plot3(x, y, z, 'b'); % Dibuja el wireframe de la cabina
    end

    % 3. Zona de transición
    posicionesTransicion = linspace(seccionTransicion(1), seccionTransicion(2), 5);
    elipsesTransicion = [];
    for pos = posicionesTransicion
        [x, y, z] = elipse3D(radioTransicion(1), radioTransicion(2), pos, elevacionTransicion, numPuntos);
        elipsesTransicion = [elipsesTransicion; x', y', z'];
        plot3(x, y, z, 'g'); % Dibuja el wireframe de la transición
    end

    % 4. Empenaje trasero elevado
    posicionesEmpenaje = linspace(seccionEmpenaje(1), seccionEmpenaje(2), 5);
    elipsesEmpenaje = [];
    for pos = posicionesEmpenaje
        [x, y, z] = elipse3D(radioEmpenaje(1), radioEmpenaje(2), pos, elevacionEmpenaje, numPuntos);
        elipsesEmpenaje = [elipsesEmpenaje; x', y', z'];
        plot3(x, y, z, 'r'); % Dibuja el wireframe de la cola
    end

    % 5. Conexión entre elipses (líneas de refuerzo estilo wireframe)
    for i = 1:numPuntos
        % Conectar líneas de la cabina a la transición
        plot3([elipsesCabina(i, 1), elipsesTransicion(i, 1)], ...
              [elipsesCabina(i, 2), elipsesTransicion(i, 2)], ...
              [elipsesCabina(i, 3), elipsesTransicion(i, 3)], 'k');

        % Conectar líneas de la transición al empenaje
        plot3([elipsesTransicion(i, 1), elipsesEmpenaje(i, 1)], ...
              [elipsesTransicion(i, 2), elipsesEmpenaje(i, 2)], ...
              [elipsesTransicion(i, 3), elipsesEmpenaje(i, 3)], 'k');
    end

    hold off;
end

function [x, y, z] = crearNariz(longitudCono, radioBase, numPuntos)
    % crearNariz - Genera las coordenadas de una nariz puntiaguda en 3D
    theta = linspace(0, 2 * pi, numPuntos); % Ángulo para la circunferencia
    xBase = radioBase(1) * cos(theta);      % Coordenadas X de la base
    yBase = radioBase(2) * sin(theta);      % Coordenadas Y de la base
    zBase = longitudCono * ones(1, numPuntos); % Coordenadas Z de la base

    % Coordenadas del vértice de la nariz (punta en el origen)
    x = [xBase, 0];
    y = [yBase, 0];
    z = [zBase, 0];
end

function [x, y, z] = elipse3D(a, b, posZ, elevacion, numPuntos)
    % elipse3D - Genera una elipse en 3D en la posición especificada
    theta = linspace(0, 2 * pi, numPuntos);
    x = a * cos(theta);                       % Coordenada X
    y = b * sin(theta) + elevacion;           % Coordenada Y (altura)
    z = posZ * ones(1, numPuntos);            % Coordenada Z (longitudinal)
end
