function [eta_0] = Xcontent(eta_0,Xc,Xm)

%% Crystal content
% This function accounts for the influence of suspensions of particles on the viscosity of viscous liquids. The model is from [5];

%% Viscosity

eta_0 = eta_0*((1-(Xc/Xm)).^-2);

end
%%
% [5] Mueller, S., Llewellin, E., and Mader, H., 2010, The rheology of suspensions of solid particles: Proceedings of the Royal Society A: Mathematical, Physical and Engineering Science, v. 466, no. 2116, p. 1201-1228.