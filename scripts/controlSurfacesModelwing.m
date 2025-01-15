function [camberAleron, AleronWingInitialLocationPorcentage, AleronWingLastLocationPorcentage] = controlSurfacesModelwing(aleronDeflection)
% La funcion pide variables como el largo de cuerda en metros, el
% porcentaje de camber en funcion de la cuerda, la localizacion del camber
% maximo y te brinda el camber nuevo alterado al deflectar el angulo del
% aleron

% Primero se obtiene la geometria necesaria con la que trabajar, al
% deflectar el angulo o "panel" de aleron, se obtiene una figura geometrica
% nueva, utilizando el teoorema de senos y cosenos se obtendran los angulos
% necesarios para determinar el camber nuevo

% El paso 1 es determinar angulos y distancias de lados del triangulo para
% oder trabajar con este. Distancia del panel 1, distancia del panel 2,
% distancia del panel del aleron en funcion de porcentaje de cuerda de
% aleron y angulo estatico en camber
% SE TRABAJA EN RADIANES

data = readtable('geometricData.txt', 'Delimiter', '\t');
AleronWingInitialLocationPorcentage = data.AleronWingInitialLocationPorcentage; % aleron location in wing in porcentage
AleronWingLastLocationPorcentage = data.AleronWingLastLocationPorcentage;
camberPorcentage = data.camberWingPorcentage; % Camber Porcentage
camberLocationPorcentage = data.camberLocationPorcentage;
aleronChordPorcentage = data.aleronChordPorcentage;
rootChord = data.rootChordWing;
tipChord = data.tipChordWing;

% ---MAC calculation for aleron chord
m = (tipChord - rootChord)/100;
cr1 = m*AleronWingInitialLocationPorcentage + rootChord;
ct1 = m*AleronWingLastLocationPorcentage + rootChord;
taper1= ct1/cr1;
chord= (2/3)*cr1*((taper1^2 + taper1 + 1)/(taper1 + 1));
camberAleron = chord*camberPorcentage/100;
aleronChord = chord*aleronChordPorcentage/100;
camberLocation = chord*camberLocationPorcentage/100;

% ---Panels analysis
p1= sqrt(camberAleron^2 + (chord/2)^2);
a= camberAleron*(aleronChord)/(camberLocation);
paleron= sqrt(a^2 + (aleronChord)^2);
p2= p1 - paleron;
AngABC= acos((chord^2 - p1^2 - p1^2)/(-2*p1*p1));
AngBCD= deg2rad(180 + aleronDeflection);

% ---Calculo de camber---
b= sqrt(p1^2 + p2^2 - 2*p1*p2*cos(AngABC));
AngBAC= asin(sin(AngABC)*p2/b);
AngBCA= asin(sin(AngABC)*p1/b);
AngACD= AngBCD - AngBCA;
newchord= sqrt(b^2 + paleron^2 - 2*b*paleron*cos(AngACD));
AngCAD= asin(sin(AngACD)*paleron/newchord);
AngBAD= AngBAC + AngCAD;
camberAleron= p1*sin(AngBAD);

end