function[rho_f] = enthalpy(p)
%% This function calculates the density of water vapour in region 3

% Enthalpy in region 3 is averaged between the p--T- dependent enthalpies
% of regions 1 and 2, using equations from [6].

% % Simplification made between 16 and 22 MPa. Next release may account for
% this.
%     sat_T = saturation_temperature(p);

H = 0.5*(H1(p)+H2(p));

if H(m) < B3a3b(p) % sub-region == 3a
    Ii = [-12, -12, -12, -12, -10, -10, -10, -8, -8, -6, -6, -6, -4, -4, -3, -2, -2, -1, -1, -1, -1, 0, 0, 1, 1, 1, 2, 2, 3, 4, 5, 8];
    Ji = [6, 8, 12, 18, 4, 7, 10, 5, 12, 3, 4, 22, 2, 3, 7, 3, 16, 0, 1, 2, 3, 0, 1, 0, 1, 2, 0, 2, 0, 2, 2, 2];
    ni = [5.29944062966028E-03, -0.170099690234461, 11.1323814312927, -2178.98123145125, -5.06061827980875E-04, 0.556495239685324, -9.43672726094016, -0.297856807561527, 93.9353943717186, 1.92944939465981E-02, 0.421740664704763, -3689141.2628233, -7.37566847600639E-03, -0.354753242424366, -1.99768169338727, 1.15456297059049, 5683.6687581596, 8.08169540124668E-03, 0.172416341519307, 1.04270175292927, -0.297691372792847, 0.560394465163593, 0.275234661176914, -0.148347894866012, -6.51142513478515E-02, -2.92468715386302, 6.64876096952665E-02, 3.52335014263844, -1.46340792313332E-02, -2.24503486668184, 1.10533464706142, -4.08757344495612E-02];
    p_star = 100; % Reference pressure in MPa.
    H_star = 2100; % Reference enthalpy in kJ kg-1.
    pi(m) = p(m)/p_star;
    oldh = H(m)/H_star;
    omega(m) = 0; % Set initial volume to zero.
    v_star = 0.0028; % Reference volume.
    for i = 1 : 32
        omega(m) = omega(m) + ni(i) * (pi(m) + 0.128) .^ Ii(i) * (oldh - 0.727) .^ Ji(i);
    end
    v3(m) = omega(m)*v_star;
else         % sub-region == 3b
    Ii = [-12, -12, -8, -8, -8, -8, -8, -8, -6, -6, -6, -6, -6, -6, -4, -4, -4, -3, -3, -2, -2, -1, -1, -1, -1, 0, 1, 1, 2, 2];
    Ji = [0, 1, 0, 1, 3, 6, 7, 8, 0, 1, 2, 5, 6, 10, 3, 6, 10, 0, 2, 1, 2, 0, 1, 4, 5, 0, 0, 1, 2, 6];
    ni = [-2.25196934336318E-09, 1.40674363313486E-08, 2.3378408528056E-06, -3.31833715229001E-05, 1.07956778514318E-03, -0.271382067378863, 1.07202262490333, -0.853821329075382, -2.15214194340526E-05, 7.6965608822273E-04, -4.31136580433864E-03, 0.453342167309331, -0.507749535873652, -100.475154528389, -0.219201924648793, -3.21087965668917, 607.567815637771, 5.57686450685932E-04, 0.18749904002955, 9.05368030448107E-03, 0.285417173048685, 3.29924030996098E-02, 0.239897419685483, 4.82754995951394, -11.8035753702231, 0.169490044091791, -1.79967222507787E-02, 3.71810116332674E-02, -5.36288335065096E-02, 1.6069710109252];
    p_star = 100; % Reference pressure in MPa.
    H_star = 2800; % Reference enthalpy in kJ kg-1.
    pi(m) = p(m)/p_star;
    oldh = H/H_star;
    omega(m) = 0; % Set initial volume to zero.
    v_star = 0.0088; % Reference volume.
    for i = 1 : 30
        omega(m) = omega(m) + ni(i) * (pi + 0.0661) .^ Ii(i) * (oldh - 0.72) .^ Ji(i);
    end
    v3(m) = omega(m)*v_star;
end
rho_f(m) = 1/v3;
end

function [H1] = H1(p)
T = 623.15; % Limit temperature in K
R = 0.461526; %kJ kg-1 K-1
Ii = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32];
Ji = [-2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41];
ni = [0.14632971213167, -0.84548187169114, -3.756360367204, 3.3855169168385, -0.95791963387872, 0.15772038513228, -0.016616417199501, 8.1214629983568E-04, 2.8319080123804E-04, -6.0706301565874E-04, -0.018990068218419, -0.032529748770505, -0.021841717175414, -5.283835796993E-05, -4.7184321073267E-04, -3.0001780793026E-04, 4.7661393906987E-05, -4.4141845330846E-06, -7.2694996297594E-16, -3.1679644845054E-05, -2.8270797985312E-06, -8.5205128120103E-10, -2.2425281908E-06, -6.5171222895601E-07, -1.4341729937924E-13, -4.0516996860117E-07, -1.2734301741641E-09, -1.7424871230634E-10, -6.8762131295531E-19, 1.4478307828521E-20, 2.6335781662795E-23, -1.1947622640071E-23, 1.8228094581404E-24, -9.3537087292458E-26];

p_star = 16.53; % Reference pressure.
pi = p/p_star; % Dimensionless pressure.
T_star = 1386; % Reference temperature.
tau = T_star/T; % Dimensionless temperature.
gamma_tau = 0;  % Set initial value to zero.
for i = 1 : 34
    gamma_tau = gamma_tau + (ni(i)*(7.1 - pi).^Ii(i) * Ji(i) * (tau - 1.222).^(Ji(i) - 1));
end
H1 = R*T*tau*gamma_tau;
end

function [H2] = H2(p)

T = B23T(p); % Limit temperature in K (boundary of Regions 2 and 3)
R = 0.461526; %kJ kg-1 K-1
Jo = [0, 1, -5, -4, -3, -2, -1, 2, 3];
no = [-9.6927686500217, 10.086655968018, -0.005608791128302, 0.071452738081455, -0.40710498223928, 1.4240819171444, -4.383951131945, -0.28408632460772, 0.021268463753307];
Ir = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24];
Jr = [0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35, 0, 11, 25, 8, 36, 13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58];
nr = [-0.17731742473213E-2, -0.17834862292358E-1, -0.45996013696365E-1, -0.57581259083432E-1, -0.05032527872793, -3.3032641670203E-05, -1.8948987516315E-04, -3.9392777243355E-03, -0.043797295650573, -2.6674547914087E-05, 2.0481737692309E-08, 4.3870667284435E-07, -3.227767723857E-05, -1.5033924542148E-03, -0.040668253562649, -7.8847309559367E-10, 1.2790717852285E-08, 4.8225372718507E-07, 2.2922076337661E-06, -1.6714766451061E-11, -2.1171472321355E-03, -23.895741934104, -5.905956432427E-18, -1.2621808899101E-06, -0.038946842435739, 1.1256211360459E-11, -8.2311340897998, 1.9809712802088E-08, 1.0406965210174E-19, -1.0234747095929E-13, -1.0018179379511E-09, -8.0882908646985E-11, 0.10693031879409, -0.33662250574171, 8.9185845355421E-25, 3.0629316876232E-13, -4.2002467698208E-06, -5.9056029685639E-26, 3.7826947613457E-06, -1.2768608934681E-15, 7.3087610595061E-29, 5.5414715350778E-17, -9.436970724121E-07];
pi = p; % p* = 1 MPa, thus pi = p.
T_star = 540; % Reference temperature
tau = T_star./T; % Dimensionless temperature

gamma_o_tau = 0; % Set initial value to zero.
gamma_r_tau = 0; % Set initial value to zero.

for i = 1:9;
    gamma_o_tau = gamma_o_tau + no(i)*Jo(i)*tau.^(Jo(i)-1);
end

for m = 1:size(p);
    for j = 1:43;
        gamma_r_tau = gamma_r_tau + nr(j)*pi(m).^Ir(j)*Jr(j)*(tau - 0.5).^(Jr(j)-1);
    end
end

H2 = R*T*tau*(gamma_o_tau + gamma_r_tau);
end

function [B23T] = B23T(p)
%% B23 function required to define boundary between regions 2 and 3 as a function of pressure.
% Based on [3]
pi = p; % p* = 1 MPa, thus pi = p.

n3 = 0.10192970039326E-2;
n4 = 0.57254459862746E3;
n5 = 0.13918839778870E2;
B23T = n4 + ((pi - n5)/n3).^0.5;
end

function [B3a3b] = B3a3b(p)
%% Boundary equation between subRegion 3a and 3b.
n1 = 0.201464004206875E4;
n2 = 0.374696550136983E1;
n3 = -0.219921901054187E-1;
n4 = 0.87513168600995E-4;
%B3a3b = NaN(size(p));

B3a3b = n1+(n2*p)+(n3*p.^2)+(n4*p.^3);
end

% function[sat_T] = saturation_temperature(p) % IAPWS 2007
% %This function calculates the saturation temperature as a function of
% %pressure of H2O. Based on [3 %2007].
% 
% % Numerical values of the coefficients of the dimensionless (temperature) saturation equation
% ni1 = 0.11670521452767E4;
% ni2 = -0.72421316703206E6;
% ni3 = -0.17073846940092E2;
% ni4 = 0.1202082470247E5;
% ni5 = -0.32325550322333E7;
% ni6 = 0.1491510861353E2;
% ni7 = -0.48232657361591E4;
% ni8 = 0.40511340542057E6;
% ni9 = -0.23855557567849;
% ni10 = 0.65017534844798E3;
% 
% zeta = p^0.25;
% E = (zeta^2)+(n3*zeta)+ n6;
% F = (n1*zeta^2) + (n4*zeta) + n7;
% G = (n2*zeta^2) + (n5*zeta) + n8;
% D = (2*G)/-F-(F^2-(4*E*G))^0.5;
% sat_T = (n10 + D -((n10+D)^2-4*(n9+n10*D)))/2; % T = saturation T/ T*, where T* = 1 K.
% end

%%
% [6] IAPWS, 2014. Revised Supplementary Release on Backward Equations for the Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam. IAPWS meeting June 2014: Moscow, Russia. 
