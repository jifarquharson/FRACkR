function [H2O] = equilibrium(T,sigma)

%% Pressure dependence of H2O solubility

% The equlibrium water content is given as a function of pressure and temperature, after Liu et al., 2005 [3].
% Pressure in MPa
p = sigma./1000000; 
% Temperature in Kelvin
T = T + 273.15; 
s1 = 354.94;
s2 = 9.623;
s3 = 1.5223;
s4 = 0.0012439;
%%
Ceq = (((s1*p.^0.5)+(s2*p)-(s3*p.^1.5))/T)+(s4*p.^1.5);
H2O = Ceq;

end

%% 
% [3] Liu, Y., Zhang, Y., and Behrens, H., 2005, Solubility of H2O in rhyolitic melts at low pressures and a new empirical model for mixed H2O–CO2 solubility in rhyolitic melts: Journal of Volcanology and Geothermal Research, v. 143, no. 1, p. 219-235.

