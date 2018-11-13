function [time, dt, T, window_max, window_min, l, nf, mean_wf, intact, spacing, w, sizew, Xc] = inputs(nt)
%% Primary input arguments
%   Prompted user inputs, which outputs parameters used in models.

prompt = 'Timescale to compute [days]? ';
time = input(prompt); dt = (time*86400)./ nt;
T = input('What is the temperature [°C]? ');

% Maximum depth of segment of interest
window_max = input('What is the lower limit of the lengthscale [m]? ');
% Minimum depth of segment of interest
window_min = input('What is the upper limit of the lengthscale [m]? ');
% Range of interest
window = window_max - window_min; l = window;
% Fracture density input
rho_frac = input('What is the fracture density [fractures/m]? ');
% Calculates total number of fractures
nf = rho_frac*window; nf = round(nf);
% Mean fracture width input
[mean_wf] = fracture_width(nf,l);
% Width of non-fractured laterial [m];
intact = window - (nf*mean_wf);
% Mean spacing in between fractures [m];
spacing = intact/nf;

% Generates fracture depths up to a user-defined density
[w] = fracture_generator(window_min, window_max,nf, mean_wf);
% Cumulative fracture width
sizew = size(w,1);

% Call crystal content input
[Xc] = Crystal();

end

