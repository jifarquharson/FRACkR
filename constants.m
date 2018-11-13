function [alpha, phi_i, rho, g, A, B, C, D, phi_c] = constants()
%% Imports constants
% Imports constants used in sintering model [1] and permeability model
% [2]. To alter constants, the values below can be changed.

%% Sintering constants
alpha = 2;
%empirical coefficient used in porosity reduction model [1]
phi_i = 0.4;
%intial porosity
rho = 2500;
%initial bulk density [kg/m^3]
g = 9.806;
%acceleration due to gravity [m^2/s]

%% Permeability-porosity relation constants
A = 7.97686606*10^-26;
% k - phi relation constant [2]
B = 8.76;
% k - phi relation constant [2]
C = 1.335903*10^-16;
% k - phi relation constant [2]
D = 1.01;
% k - phi relation constant [2]
phi_c = 0.155;
% changepoint porosity [2]
end

%%
% [1] Russell, J. K., & Quane, S. L. (2005). Rheology of welding: inversion of field constraints. Journal of Volcanology and Geothermal Research, 142(1), 173-191.
%%
% [2] Heap, M. J., Farquharson, J. I., Wadsworth, F. B., Kolzenburg, S., & Russell, J. K. (2015). Timescales for permeability reduction and strength recovery in densifying magma. Earth and Planetary Science Letters, 429, 223-233.
