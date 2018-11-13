function[rho_f]=density(sigma,T, region)
%% This function calculates the density of H2O as a function of pressure, temperature and p--T region, after [6] and [8].

R = 0.461526; % Specific gas constant [kJ kg^-1 K^-1]
p = sigma./1000000; % Pressure in MPa
T = T + 273.15;      % Temperature in Kelvin
rho_f = NaN(size(p)); % Pre-allocate empty matrix

for m = 1:size(p)
    if region == 1
        I = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32];
        J = [-2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41];
        n = [0.14632971213167, -0.84548187169114, -3.756360367204, 3.3855169168385, -0.95791963387872, 0.15772038513228, -0.016616417199501, 8.1214629983568E-04, 2.8319080123804E-04, -6.0706301565874E-04, -0.018990068218419, -0.032529748770505, -0.021841717175414, -5.283835796993E-05, -4.7184321073267E-04, -3.0001780793026E-04, 4.7661393906987E-05, -4.4141845330846E-06, -7.2694996297594E-16, -3.1679644845054E-05, -2.8270797985312E-06, -8.5205128120103E-10, -2.2425281908E-06, -6.5171222895601E-07, -1.4341729937924E-13, -4.0516996860117E-07, -1.2734301741641E-09, -1.7424871230634E-10, -6.8762131295531E-19, 1.4478307828521E-20, 2.6335781662795E-23, -1.1947622640071E-23, 1.8228094581404E-24, -9.3537087292458E-26];
        p_star = 16.53; % Reference pressure in MPa
        pi(m) = p(m)/p_star; % dimensionless pressure [p/p*]
        T_star = 1386; % Reference temperature. 
        tau = T_star/T; % inverse dimensionless temperature [T*/T] [Region 1]
        oldg_pi = 0; % Initial condition for specific Gibbs free energy (dimensionless)
        for i = 1:34
            oldg_pi = oldg_pi - n(i) * I(i) * (7.1 - pi) .^ (I(i) - 1) * (tau - 1.222) .^ J(i);
        end
        v1(m) = R * T / p(m) * pi(m) * oldg_pi / 1000; % Specific volume
        rho_f(m) = 1 / v1(m); % Gas density for region 1.
        
    elseif region == 2
        Ir = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24];
        Jr = [0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35, 0, 11, 25, 8, 36, 13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58];
        nr = [-1.7731742473213E-03, -0.017834862292358, -0.045996013696365, -0.057581259083432, -0.05032527872793, -3.3032641670203E-05, -1.8948987516315E-04, -3.9392777243355E-03, -0.043797295650573, -2.6674547914087E-05, 2.0481737692309E-08, 4.3870667284435E-07, -3.227767723857E-05, -1.5033924542148E-03, -0.040668253562649, -7.8847309559367E-10, 1.2790717852285E-08, 4.8225372718507E-07, 2.2922076337661E-06, -1.6714766451061E-11, -2.1171472321355E-03, -23.895741934104, -5.905956432427E-18, -1.2621808899101E-06, -0.038946842435739, 1.1256211360459E-11, -8.2311340897998, 1.9809712802088E-08, 1.0406965210174E-19, -1.0234747095929E-13, -1.0018179379511E-09, -8.0882908646985E-11, 0.10693031879409, -0.33662250574171, 8.9185845355421E-25, 3.0629316876232E-13, -4.2002467698208E-06, -5.9056029685639E-26, 3.7826947613457E-06, -1.2768608934681E-15, 7.3087610595061E-29, 5.5414715350778E-17, -9.436970724121E-07];
        pi(m) = p(m); % dimensionless pressure (p/p*; p* = 1) [Region 2]
        T_star = 540; % Reference temperature.
        tau = T_star/T; % inverse dimensionless temperature [T*/T] [Region 2]
        oldg_o_pi = 1/pi(m); 	% Initial condition for specific Gibbs free energy (dimensionless, Ideal Gas part)
        oldg_r_pi = 0;	% Initial condition for specific Gibbs free energy (dimensionless, residual part)
        for i = 1:43 
            oldg_r_pi = oldg_r_pi + nr(i)*Ir(i)*pi(m).^(Ir(i) - 1)*(tau - 0.5).^Jr(i);
        end
        v2(m) = R*T/p(m)*pi(m)*(oldg_o_pi + oldg_r_pi)./1000;
        rho_f(m) = 1 / v2(m); % Gas density for region 2.
        
    elseif region == 3
        rho_f(m) = enthalpy(p); % Calls enthalpy function, which approximates fluid density as a function of enthalpy.

    elseif region == 4
        rho_f(m) = NaN;

    elseif region == 5
        Ir = [1, 1, 1, 2, 3];
        Jr = [0, 1, 3, 9, 3];
        nr = [-1.2563183589592E-04, 2.1774678714571E-03, -0.004594282089991, -3.9724828359569E-06, 1.2919228289784E-07];
        T_star = 1000; % Reference temperature.
        tau =  T_star/T; % inverse dimensionless temperature [T*/T] [Region 5]
        pi(m) = p(m); % dimensionless pressure (p/p*; p* = 1) [Region 5]
        oldg_o_pi = 1/pi(m);
        oldg_r_pi = 0;
        for i = 1:5
            oldg_r_pi = oldg_r_pi + nr(i)*Ir(i)*pi(m).^(Ir(i) - 1)*tau .^ Jr(i);
        end
        v5(m) = R*T/p(m)*pi(m)*(oldg_o_pi + oldg_r_pi)/1000;
        rho_f(m) = 1/v5(m); % Gas density for region 5.
    end
rho_f = rho_f';
end

%%
% [6] IAPWS, 2014. Revised Supplementary Release on Backward Equations for the Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam. IAPWS meeting June 2014: Moscow, Russia. 

% [8] IAPWS, 2012. Revised Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam. IAPWS meeting August 2007: Lucerne, Switzerland.