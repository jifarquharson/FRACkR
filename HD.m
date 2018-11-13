function [eta_0] = HD(T,H2O,Xc)

%% This viscosity model for calc-alkaline rhyolites is from Hess and Dingwell, 1996 [4]

%% Constants

a1 = -3.545;
a2 = 0.833;
b1 = 9601;
b2 = -2368;
c1 = 195.7;
c2 = 32.25;

%%
a = a1 + (a2*log(H2O));
b = b1 + (b2*log(H2O));
c = c1 + (c2*log(H2O));


%% Viscosity calculation

l_eta0 = a+(b./((T+273.15)-c));
eta_0 = 10.^l_eta0;

rp = 1;
b = 1.08;

Xm_x = 0.656;
Xm = Xm_x*exp(-((log10(rp))^2)/(2*b*b));

%% Call particle suspension function
[eta_0] = Xcontent(eta_0,Xc,Xm);

end
%%
% [4] Hess, K. U., and Dingwell, D. B., 1996, Viscosities of hydrous leucogranitic melts: A non-Arrhenian model: American Mineralogist, v. 81, no. 9-10, p. 1297-1300.