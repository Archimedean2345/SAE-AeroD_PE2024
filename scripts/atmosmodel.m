function [T, P, rho, mu, SpeedSound] = atmosmodel(altitude,localTemperature,localAltitude)
    % Constants
    g0= 9.80665;          % Sea level gravitational acceleration (m/s^2)
    R= 287.05;            % Specific gas constant for dry air (J/(kg*K))
    P0= 101325;           % Sea level standard pressure (Pa)
    B= 0.0065;            % Temperature lapse rate (K/m) in the troposphere
    alt_tropopause= 11000; % Altitude of the tropopause (m)
    % Temperature correction
    T0= (localTemperature+ 273.15) + B*localAltitude;
     % Sutherland's law constants for air
    mu0 = 1.716e-5;        % Reference viscosity at T0 (Pa.s)
    C = 110.4;             % Sutherland's constant (K)
    gamma= 1.4;             % heat capacity ratio for air

    % Troposphere (0 to 11 km)
    if altitude <= alt_tropopause
        T= T0 - B*altitude;
        P= P0*(T/T0)^(g0/(R*B));
        rho= P/(R*T);
        mu= mu0*((288.15+ C)/(T+ C))*(T/288.15)^(3/2);
        SpeedSound= sqrt(gamma*R*T);
    % Stratosphere (11 km to 20 km)
    elseif altitude <= 20000
        T= 216.65; % Constant temperature in the stratosphere (K)
        P= P0*(216.65/T0)^(g0/(R*B))*exp(-g0*(altitude-alt_tropopause)/(R*T));
        rho= P/(R*T);
        mu= mu0*((288.15+ C)/(T+ C))*(T/288.15)^(3/2);
        SpeedSound= sqrt(gamma*R*T);
    else
        error('Altitude out of range. The model is valid from 0 to 20 km.');
    end
end