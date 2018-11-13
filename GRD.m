function [eta_0] = GRD(T,H2O,Xc)

%% This viscosity model for volatile-bearing melts is from Giordano et al. (2008) and should be cited as such. 
%Compositional inputs required are the oxide fractions of SiO2, TiO2, Al2O3, FeO, MnO, MgO, CaO, Na2O, K2O, P2O5, and F2O-1, in weigth percent (wt.%). 

%% Inputs

SiO2 = input('SiO2 content [wt.%]? ');
TiO2 = input('TiO2 content [wt.%]? ');
Al2O3 = input('Al2O3 content [wt.%]? ');
FeO = input('FeO content [wt.%]? ');
MnO = input('MnO content [wt.%]? ');
MgO = input('MgO content [wt.%]? ');
CaO = input('CaO content [wt.%]? ');
Na2O = input('Na2O content [wt.%]? ');
K2O = input('K2O content [wt.%]? ');
P2O5 = input('P2O5 content [wt.%]? ');
%H2O = input('H2O content [wt.%]? ');
F2O = input('F2O-1 content [wt.%]? ');

Oxides = [SiO2;TiO2;Al2O3;FeO;MnO;MgO;CaO;Na2O;K2O;P2O5;H2O;F2O];
Oxides = Oxides./(sum(Oxides))*100;
% Normalised wt.% oxides

%% Constants

Molwt = [60.085;79.88;101.96;71.85;70.94;40.30;56.08;61.98;94.20;141.94;18.001;37.9968]; % Molecular weight
gfw = 61.0876;
% Gram fresh weight
OB = gfw.*Oxides./Molwt;
% Mol.% oxide basis

A = -4.55;

%%
B1 = 159.60*(OB(1)+OB(2));
B2 = -173.3*(OB(3));
B3 = 72.10*((OB(4))+(OB(5))+(OB(10)));
B4 = 75.70*OB(6);
B5 = -39.00*OB(7);
B6 = -84.10*(OB(8)+OB(11)+OB(12));
B7 = 141.50*(OB(11)+OB(12)+(log(1+(OB(11)))));
B11 = -2.43*(OB(1)+OB(2))*(OB(4)+OB(5)+OB(6));
B12 = -0.91*(OB(1)+OB(2)+OB(3)+OB(10))*(OB(8)+OB(9)+OB(11));
B13 = 17.60*(OB(3))*(OB(8)+OB(9));

B =(B1+B2+B3+B4+B5+B6+B7+B11+B12+B13);

%%
C1 = 2.75*OB(1);
C2 = 15.70*(OB(2)+OB(3));
C3 = 8.30*(OB(4)+OB(5)+OB(6));
C4 = 10.20*OB(7);
C5 = -12.30*(OB(8)+OB(9));
C6 = -99.50*(log(1+OB(11)+OB(12)));
C11 = 0.30*(OB(3)+OB(4)+OB(5)+OB(6)+OB(7)+OB(10))*(OB(8)+OB(9)+OB(11)+OB(12));

C = (C1+C2+C3+C4+C5+C6+C11);

%% Viscosity calculation

log_eta0 = A+(B/((T+273.15)-C));
eta_0 = 10.^log_eta0;

rp = 1;
b = 1.08;

Xm_x = 0.656;
Xm = Xm_x*exp(-((log10(rp))^2)/(2*b*b));

%% Call particle suspension function
[eta_0] = Xcontent(eta_0,Xc,Xm);
end

% Giordano D, Russell JK, and Dingwell DB (2008) Viscosity of Magmatic Liquids: A
% Model. Earth and Planetary Science Letters.