function [eta_0, sigma] = viscosity(T, Xc, rho, g, w, sizew)
%% Incorporates water content into viscosity
% Allows the  user to choose between equilibrium water content or a defined
% value.
sigma = rho*g*w;
    % Calculates magmastatic (or lithostatic) stress driving densification
    
    choice = questdlg('Would you like to calculate equilibrium H2O content or input value?', ...
        'Continue', ...
        'Manual input','Equilibrium','Manual input');
    % Handle response
    switch choice
        case 'Equilibrium'
            
            for i = 1:sizew
    
    [H2O] = equilibrium(T,sigma);
    % Runs solubility function
    
            end
         case 'Manual input'
            [H2O] = input('What is the dissolved water content? ');
            
    end
    % By default, this model employs the viscosity model of Hess and
    % Dingwell (1996). This can be substituted for the model of Giordano et
    % al. (2008), by commenting and uncommenting the models below as
    % appropriate. Note that the GRD model requires the oxide fractions of
    % the melt to be known.

[eta_0] = HD(T,H2O, Xc); % Viscosity model of Hess and Dingwell
% [eta_0] = GRD(T,H2O, Xc); % Viscosity model of Giordano et al.
    % Runs viscosity function
end

