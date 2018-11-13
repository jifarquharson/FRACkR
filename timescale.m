function [delta, nu] = timescale(T, w, A, B, phi_i, eta_0, sigma)
%% Delta calculator
% This function computes the value of delta, ratio between fracture width and compaction lengthscale

[nu] = pore_fluid(sigma, T);
% Calculates viscosity of pore fluid at given temperature

%% Constants

kfj = A*((100*phi_i).^B); % Fracture permeability at intial porosity

%% Delta calculation

delta = w./((kfj*(eta_0/nu)).^0.5);

end