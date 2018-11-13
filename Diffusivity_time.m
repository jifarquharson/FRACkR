function [lambda_D] = Diffusivity_time(H2O,T,a)
% Calculates the timescale of molecular diffusion using the diffusivity of
% a rhyolitic melt (after [9]) when the characterisitc particle radius a is
% known.
%   
C0 = 1; % Reference water content in wt.%.
T = T + 273.15; % Temperature in Kelvin
DH2O = (H2O./C0)*exp(-16.83 - (10992/T)); % Diffusivity after [9]

lambda_D = (a^2)./DH2O; % Timescale of diffusion in seconds.

end

%%
% [9] Zhang, Y., Stolper, E.M. and Wasserburg, G.J., 1991. Diffusion of water in rhyolitic glasses. Geochimica et Cosmochimica Acta, 55(2), pp.441-456.