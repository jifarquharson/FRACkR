function [nu] = pore_fluid(sigma, T)
%% Calculates the viscosity of water vapour under a given p and T.
% Employs the constants and equations of IAPWS [7] and [8].

T = T + 273.15;         % Temperature in Kelvin      
p = sigma./1000000;     % Pressure in MPa

T_star = 647.096;       % Reference temperature in Kelvin after [7]
rho_star = 322.0;       % Reference density [kg m^-3] after [7]
p_star = 22.064;        % Reference pressure [MPa] after [7]
nu_star = 1.00*10^-6;   % Reference viscosity [Pa.s] after [7]

rho_f = density(sigma,T, region(sigma,T));   % Calls density function, and region function
theta = T/T_star;       	 % Dimensionless temperature
rho_var = rho_f/rho_star; 	 % Dimensionless density
pi = p/p_star;       		 % Dimensionless pressure

%% The contribution to viscosity due to finite density [nu1]

H0 = [0.520094, 0.0850895, -1.08374, -0.289555, 0, 0];
H1 = [0.222531, 0.999115, 1.88797, 1.26613, 0, 0.120573];
H2 = [-0.281378, -0.906851, -0.772479, -0.489837, -0.257040, 0];
H3 = [0.161913, 0.257399, 0, 0, 0, 0];
H4 = [-0.0325372, 0, 0, 0.0698452, 0, 0];
H5 = [0, 0, 0, 0, 0.00872102, 0];
H6 = [0, 0, 0, -0.00435673, 0, -0.000593264];

Sum = 0;
for i = 0 : 5
Sum = Sum + H0(i+1)*(1/theta - 1).^i + H1(i+1)*(1 / theta - 1).^i*(rho_var - 1).^1 + H2(i+1)*(1/theta - 1).^i*(rho_var - 1).^2 + H3(i+1)*(1/theta - 1).^i*(rho_var - 1).^3 + H4(i+1)*(1/theta - 1).^i*(rho_var - 1).^4 + H5(i+1)*(1/theta - 1).^i*(rho_var - 1).^5 + H6(i+1)*(1/theta - 1).^i*(rho_var - 1).^6;
end

nu1 = exp(rho_var.*Sum);

%% Viscosity in the dilute-gas limit [nu0]

Hi0 = 1.67752;
Hi1 = 2.20462;
Hi2 = 0.6366564;
Hi3 = -0.241605;

nu0 = (100.*sqrt(theta))/((Hi0)+(Hi1/theta)+(Hi2/theta^2)+(Hi3/theta^3));

%% Pore fluid (water) viscosity as a function of temperature and pressure[Pa.s]

Omega = nu1*nu0;
nu = Omega*nu_star;
end
%% 
% [7] IAPWS, 2008. Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary Water Substance. IAWPS meeting September 2008: Berlin, Germany.
% [8] IAPWS, 2012. Revised Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam. IAPWS meeting August 2007: Lucerne, Switzerland.
