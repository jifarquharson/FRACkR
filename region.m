function [region] = region(sigma,T)
%% Region selection in p--T space
%Based on IAPWS [6] and [8]. 

T = T + 273.15;      % Temperature in Kelvin
p = sigma./1000000; % Pressure in MPa

%   +-------+---+-----+
%   |       |   |     |
%   |       |   |     |
%   |       |   |     |
%   |       |   |     +----------------+
% p |   1   | 3 |  2  |                |
%   |       |  /      |                |
%   |       |_/       |        5       |
%   |      /4         |                |
%   +----=-----------------------------+
%                   T

region = NaN(size(p)); % Pre-allocate empty matrix

for i = 1: size(p)
% Selects region (1 - 5, as above) as a function of pressure p (MPa) and temperature T (K)
if T > 1073.15 && T < 2273.15 &&  p(i) <= 50 && p(i) > 0.000611 
    region(i) = 5; % p - T domain = Region 5
elseif T <= 1073.15 && T > 273.15 && p(i) <= 100 && p(i) > 0.000611 
    if T > 623.15 
        if p(i) > B23(T) % p must be above the B23 line as a function of temperature
            region(i) = 3; % p - T domain of Region 3    
            if T < 647.096 
                sat_p = saturation_pressure(T);
                if abs(p(i) - sat_p) < 0.00001 % Only the small region around the saturation pressure is defined here as region 4.
                    region(i) = 4;
                end
            end
        else
            region(i) = 2; % all other p - T relevant conditions above 623.15 K are defined as Region 2
        end
    else    % defines conditions in low-temperature domain
        sat_p = saturation_pressure(T);
        if abs(p(i) - sat_p) < 0.00001 % Region 4 defined as above
            region(i) = 4;
        elseif p(i) > sat_p
            region(i) = 1; % Region 1 defined as the region above the saturation pressure
        else
            region(i) = 2; % extension of Region 2 below the saturation pressure.
        end
    end
else
    region(i) = 0; % Imposed p - T conditions outwith valid p - T domain.
end
region = region';
end
end

function [B23] = B23(T)
%% B23 function required to define boundary between regions 2 and 3 as a function of temperature.
% Based on [3]
T = T + 273.15; % Temperature in Kelvin

n1 = 0.34805185628969E3;
n2 = -0.11671859879975E1;
n3 = 0.10192970039326E-2; 
B23 = n1 + (n2*T) + (n3*T^2);
end

function [sat_p] = saturation_pressure(T)
%This function calculates the saturation pressure as a function of
%temperature of H2O. Based on [8].

T = T + 273.15; % Temperature in Kelvin

% Numerical values of the coefficients of the dimensionless (pressure) saturation equation
n1 = 0.11670521452767E4;
n2 = -0.72421316703206E6;
n3 = -0.17073846940092E2;
n4 = 0.1202082470247E5;
n5 = -0.32325550322333E7;
n6 = 0.1491510861353E2;
n7 = -0.48232657361591E4;
n8 = 0.40511340542057E6;
n9 = -0.23855557567849;
n10 = 0.65017534844798E3;

theta = T + n9/(T - n10);

A = theta^2 + (n1 * theta) + n2;
B = (n3*theta^2) +( n4 * theta) + n5;
C = (n6 * theta^2) + (n7 * theta) + n8;

% Saturation pressure
sat_p = (2*C/(-B + (B^2-4*A*C)^0.5))^4;

end
%% 
% [6] IAPWS, 2014. Revised Supplementary Release on Backward Equations for the Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam. IAPWS meeting June 2014: Moscow, Russia. 
% [8] IAPWS, 2012. Revised Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam. IAPWS meeting August 2007: Lucerne, Switzerland.