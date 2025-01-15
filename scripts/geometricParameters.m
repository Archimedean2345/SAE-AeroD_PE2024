% ---- Geometric Parameters ---
function [AR, AoALocal, LocalChord, SurfaceLocal, SurfaceTotal, TipChord, RootChord, Span, Segments, camberPorcentage] = geometricParameters()
    % --- Read Geometric Parameters from File ---
    data = readtable('geometricData.txt', 'Delimiter', '\t');
    AR = data.ARwing; % Aspect Ratio
    SurfaceTotal = data.surfaceWing; % Total Surface Area of the Wing
    TipChord = data.tipChordWing; % Tip Chord Length
    RootChord = data.rootChordWing; % Root Chord Length
    Span = data.spanWing; % Wing Span
    camberPorcentage = data.camberWingPorcentage; % Camber Percentage
    IncidenceAngle = data.incidenceAngleWing; % Wing Incidence Angle
    TwistAngle = data.twistAngleWing; % Wing Twist Angle

    % Additional Parameters
    Segments = 50; % Number of segments, can be modified as needed

    % Initialize arrays to store results for each segment
    span_positions = linspace(0, Span / 2, Segments + 1); % Position along the semi-span (includes an additional segment)
    AoALocal = zeros(1, Segments + 1); % +1 for the last segment
    LocalChord = zeros(1, Segments + 1);
    SurfaceLocal = zeros(1, Segments);

    % Calculate local chord and local angle of attack for each position
    for i = 1:(Segments + 1)
        spanPosition = span_positions(i); % Position of the segment along the semi-span
        AoALocal(i) = (TwistAngle / (Span / 2)) * spanPosition + IncidenceAngle;
        LocalChord(i) = ((TipChord - RootChord) / (Span / 2)) * spanPosition + RootChord;
    end

    % Calculate the local surface area using the trapezoidal rule
    for i = 1:Segments
        h = span_positions(i + 1) - span_positions(i); % Segment length
        SurfaceLocal(i) = ((LocalChord(i) + LocalChord(i + 1)) / 2) * h;
    end

    % Update positions for midpoint calculations of AoA and chord length
    span_positions = linspace(Span / 2 / Segments - (Span / 2 / Segments / 2), (Span / 2) - (Span / 2 / Segments / 2), Segments);
    AoALocal = zeros(1, Segments);
    LocalChord = zeros(1, Segments);

    % Corrected local chord and AoA for each midpoint
    for i = 1:Segments
        spanPosition = span_positions(i); % Position of the segment along the semi-span
        AoALocal(i) = (TwistAngle / (Span / 2)) * spanPosition + IncidenceAngle;
        LocalChord(i) = ((TipChord - RootChord) / (Span / 2)) * spanPosition + RootChord;
    end
end
