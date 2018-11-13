%%

% The following script allows the calculation of the equivalent permeability
% ke of the wall of a volcanic conduit containing multiple fractures. The
% time-dependence of ke is captured by considering depth- and temperature-
% dependent fracture permeability reduction by viscous sintering.


clearvars; close all; clc

%% Create file

%Filename for an the output files, which save equivalent permeability
%evolution with time.
xlsFilename = input('Enter the name for the xls file: ', 's');

%% Model setup

% Intial time
t0 = 0;
% Number of timesteps to compute
nt = 100;

% Call constants function
[alpha, phi_i, rho, g, A, B, C, D, phi_c] = constants();

% Call user inputs function
[time, dt, T, window_max, window_min, l, nf, mean_wf, intact, spacing, w, sizew, Xc] = inputs(nt);

% Host rock permeabilities [m^2]
n = (11:1:22);k0 = 10.^(-n); n = 12;

% Call viscosity functions
[eta_0, sigma] = viscosity(T, Xc, rho, g, w, sizew);

% Flag if initial viscosity is greater than viscosity at Tg
if eta_0 > 1e12;
    disp('Fracture material may be below Tg. This may result in non-physical values. ');
    choice = questdlg('Fracture material may be below Tg. Would you like to continue?', ...
        'Continue', ...
        'Yes','No','No');
    % Handle response
    switch choice
        case 'Yes'
            flag = 2;
        case 'No'
            flag = 1;
    end
end
%% Preallocate matrices

% Creates a matrix to store porosity results
P = NaN(nt,nf);
% Creates a matrix to store time
tt = zeros(nt,nf);
% Creates a matrix with host rock permeabilities
kk = repmat(k0,nt,1);
% Creates a matrix to store permeability data
KK = NaN(nt,nf);
% Creates a matrix to store equivalent permeability values
ke = NaN(size(kk));

%% Check flags
% This section flags if dimensionless numbers exceed imposed thresholds

% Determines delta, the ratio  between fracture width and compaction lengthscale
[delta, nu] = timescale(T, mean_wf, A, B, phi_i, eta_0, sigma);
% Determines the Darcy compaction number,
[Dac, kcr, flag] = Darcy_compaction(eta_0, A, B, phi_i, alpha, T, mean_wf, sigma);


if any(delta) > 1e0;
    disp('One or more fractures may be too wide to sinter, or may necessitate more complex sintering criteria. This may result in non-physical values. ');
    
    choice = questdlg('One or more fractures may be too wide to sinter. Would you like to continue?', ...
        'Continue', ...
        'Yes','No','No');
    % Handle response
    switch choice
        case 'Yes'
            flag = 2;
        case 'No'
            flag = 1;
    end
end

%% Porosity - permeability evolution over time

% Set time to start time
t = t0;
timespan = 10.^(linspace(0,log10(time*86400),nt)).';

for m = 1:nt
    
    % porosity starts at intial value of 0.4, a typical porosity for polydisperse granular materials close to their maximum packing.
    phi = phi_i;
    % time increases by dt increment each iteration
    phi0 = phi;
    t = timespan(m);
    
    for i = 2:nt
        
        % Porosity evolution with time, after the model of [2]
        beta = log(((alpha.*sigma*t)./(eta_0.*(1 - phi_i))) + exp(-alpha*phi_i/(1-phi_i)));
        beta(beta>0) = 0;
        phi = beta./(beta-alpha);
        
    end
    P(m,:)= phi; P(P<0) = 0;
    
    % Permeability-porosity model described by two-slope model, after [1].
    if phi>=phi_c
        k = A*((100*phi).^B);
    else k = C*((100*phi).^D);
    end
    
    % Permeability and time matrices update with each iteration
    k(k < 0) = 0;
    KK(m,:) = k;
    tt(m,:) = t;
end

%% Quit script if user selects not to continue

if flag == 1
    error('Model quit by user. ')
end
%%
if flag == 0||2
    
    %% Equivalent permeability calculation
    
    % Intact width
    wi =(l - nf*mean_wf);wx = nf*mean_wf;
    % Fracture width times fracture permeability
    wfkf = wx.*KK; wfkf = sum(wfkf,2); wy = repmat(wfkf,1,n);
    % Intact width for each host rock permeability
    wix = repmat(wi,nt,n);
    % Equivalent permeability for host rock permeability
    ke = ((wix.*kk)+(wy))/l;
    
    %% Assess resolution
    
    if ke(nt,12)/k0(1,12)>1;
        choice = questdlg('Host material may not have returned to pre-fracture permeability. Consider re-running model using longer "Timescale to compute".', ...
            '', ...
            'Continue','Continue');
        % Handle response
        switch choice
            case 'Continue'
        end
    end
    if ke(nt*0.9,12)/k0(1,12)<1;
        choice = questdlg('For better resolution, consider re-running model using shorter "Timescale to compute".', ...
            '', ...
            'Continue','Continue');
        % Handle response
        switch choice
            case 'Continue'
        end
    end
    
    %% Figure plotting
    
    [hFig] = figures(tt, P, time, KK, t, n, ke, w, kk, window_max);
    %%
    X = ['Data are saved as ', xlsFilename, '.xls.']; msgbox(X);
    
    %% Output file
    
    % Saves user inputs on the first sheet, time and equivalent permeability
    % (in accordance with the initial host rock permeabilities) on the second,
    % and the randomly generated fracture positions in the third.
    rho_frac = nf/l;
    User_inputs = {'Max. depth [m]', window_max; 'Min. depth [m]', window_min; 'Fracture density [fractures/m]', rho_frac; 'Fracture width [m]', mean_wf; 'Temperature [C]', T; 'Crystal content', Xc};
    header1 = {'User inputs'};
    header2 = {'Time [s]', 'Equivalent permeability for host permeability of 10^-11[m2]', '10^-12[m]','10^-13[m]','10^-14[m]','10^-15[m]','10^-16[m]','10^-17[m]','10^-18[m]','10^-19[m]','10^-20[m]','10^-21[m]','10^-22[m]',};
    header3 = {'Fracture depths [m]'};
    output = [tt(:,1),ke];
    xlswrite(xlsFilename, header1);
    xlswrite(xlsFilename, User_inputs,'Sheet1','A2');
    xlswrite(xlsFilename, header2,'Sheet2');
    xlswrite(xlsFilename, output,'Sheet2','A2');
    xlswrite(xlsFilename, header3,'Sheet3');
    xlswrite(xlsFilename, w,'Sheet3','A2');
end

%% References
%%
% [2] Russell, J. K., & Quane, S. L. (2005). Rheology of welding: inversion of field constraints. Journal of Volcanology and Geothermal Research, 142(1), 173-191.
%%
% [1] Heap, M. J., Farquharson, J. I., Wadsworth, F. B., Kolzenburg, S., & Russell, J. K. (2015). Timescales for permeability reduction and strength recovery in densifying magma. Earth and Planetary Science Letters, 429, 223-233.
