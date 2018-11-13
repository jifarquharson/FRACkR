function [Dac, kcr, flag] = Darcy_compaction(eta_0, A, B, phi_i, alpha, T, mean_wf, sigma)
%% Critical permeability threshold
% This function computes the timescale ratio between compaction and Darcian flow.

flag = 0;

% Fracture permeability at intial porosity
kfj = A*((100*phi_i).^B);
% Call pore fluid viscosity function
[nu] = pore_fluid(sigma, T);
% Darcy compaction number
Dac = eta_0*kfj./(alpha*nu*(mean_wf.^2));
%%
% Critical permeability
kcr = (alpha*nu*(mean_wf/2).^2)./eta_0;
%%
if Dac < 1e0;
    disp('Pore pressure may not be in equilibrium. This may result in non-physical values.');
    
    choice = questdlg('Would you like to continue?', ...
        '', ...
        'Yes','No','No');
    % Handle response
    switch choice
        case 'Yes'
            flag = 2;
        case 'No'
            flag = 1;
    end
end
end