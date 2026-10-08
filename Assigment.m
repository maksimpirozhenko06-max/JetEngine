clear all;close all;clc;
warning off
%% To make sure that matlab will find the functions. You must change it to your situation 
relativepath_to_generalfolder='General'; % relative reference to General folder (assumes the folder is in you working folder)
addpath(relativepath_to_generalfolder); 
%% Load Nasadatabase
TdataBase=fullfile('General','NasaThermalDatabase');
load(TdataBase);
%% Nasa polynomials are loaded and globals are set. 
%% values should not be changed. These are used by all Nasa Functions. 
global Runiv Pref
Runiv=8.314472;
Pref=1.01235e5; % Reference pressure, 1 atm!
Tref=298.15;    % Reference Temperature
%% Some convenient units
kJ=1e3;kmol=1e3;dm=0.1;bara=1e5;kPa = 1000;kN=1000;kg=1;s=1;
%% Given conditions. 
%  For the final assignment take the ones from the specific case you are supposed to do.                  
v1=200;Tamb=300;P3overP2=7;Pamb=100*kPa;mfurate=0.58*kg/s;AF=204.42;             % These are the ones from the book
cFuel='H2';                                                           % Pick Gasoline as the fuel (other choices check Sp.Name)
%% Select species for the case at hand
iSp = myfind({Sp.Name},{cFuel,'O2','CO2','H2O','N2'});                      % Find indexes of these species
SpS=Sp(iSp);                                                                % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSp = length(SpS);
Mi = [SpS.Mass];
%% Air composition
Xair = [0 0.21 0 0 0.79];                                                   % Order is important. Note that these are molefractions
MAir = Xair*Mi';                                                            % Row times Column = inner product 
Yair = Xair.*Mi/MAir;                                                       % Vector. times vector is Matlab's way of making an elementwise multiplication
%% Fuel composition
Yfuel = [1 0 0 0 0];                                                        % Only fuel
%% Range of enthalpies/thermal part of entropy of species
TR = [200:1:3000];NTR=length(TR);
for i=1:NSp                                                                 % Compute properties for all species for temperature range TR 
    hia(:,i) = HNasa(TR,SpS(i));                                            % hia is a NTR by 5 matrix
    sia(:,i) = SNasa(TR,SpS(i));                                            % sia is a NTR by 5 matrix
end
hair_a= Yair*hia';                                                          % Matlab 'inner product': 1x5 times 5xNTR matrix muliplication, 1xNTR resulT -> enthalpy of air for range of T 
sair_a= Yair*sia';                                                          % same but this thermal part of entropy of air for range of T
% whos hia sia hair_a sair_a                                                  % Shows dimensions of arrays on commandline
%% Two methods are presented to 'solve' the conservation equations for the Diffusor
%-------------------------------------------------------------------------
% ----> This part shows the interpolation method
% Bisection is in the next 'cell'
%-------------------------------------------------------------------------
% [1-2] Diffusor :: Example approach using INTERPOLATION
cMethod = 'Interpolation Method';
sPart = 'Diffusor';
T1 = Tamb;
P1 = Pamb;
Rg = Runiv/MAir;
for i=1:NSp
    hi(i)    = HNasa(T1,SpS(i));
end
h1 = Yair*hi';
v2 = 0;
h2 = h1+0.5*v1^2-0.5*v2^2;                                                  % Enhalpy at stage: h2 > h1 due to kinetic energy
T2 = interp1(hair_a,TR,h2);                                                 % Interpolate h2 on h2air_a to approximate T2. Pretty accurate
for i=1:NSp
    hi2(i)    = HNasa(T2,SpS(i));
    si1(i)    = SNasa(T1,SpS(i));
    si2(i)    = SNasa(T2,SpS(i));
end
h2check = Yair*hi2';                                                        % Single value (1x5 times 5x1). Why do I do compute this h2check value? Any ideas?
s1thermal = Yair*si1';
s2thermal = Yair*si2';
lnPr = (s2thermal-s1thermal)/Rg;                                            % ln(P2/P1) = (s2-s1)/Rg , see lecture (s2 are only the temperature integral part of th eentropy)
Pr = exp(lnPr);
P2 = P1*Pr;
S1  = s1thermal - Rg*log(P1/Pref);                                          % Total specific entropy
S2  = s2thermal - Rg*log(P2/Pref);
% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,1,2);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T1,T2);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P1/kPa,P2/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v1,v2);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h1/kJ,h2/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S1/kJ,S2/kJ);
T2int = T2;

%% Two methods are presented to 'solve' the conservation equations for the Diffusor
%-------------------------------------------------------------------------
% ----> This part shows the Bisection method
%-------------------------------------------------------------------------
% [1-2] Diffusor :: Example approach using bisection (https://en.wikipedia.org/wiki/Bisection_method)
cMethod = 'Bisection Method';
sPart = 'Diffusor';
T1 = Tamb;
P1 = Pamb;
Rg = Runiv/MAir;
for i=1:NSp
    hi(i)    = HNasa(T1,SpS(i));
end
h1 = Yair*hi';
v2 = 0;
h2 = h1+0.5*v1^2-0.5*v2^2;                                                  % Enhalpy at stage: h2 > h1 due to kinetic energy
TL = T1;
TH = 1000;                                                                  % A guess for the TH (must be too high)
iter = 0;
while abs(TH-TL) > 0.01
    iter = iter+1;
    Ti = (TL+TH)/2;
    for i=1:NSp
        hi2(i)    = HNasa(Ti,SpS(i));
    end
    h2i = Yair*hi2';                                                        % Single value (1x5 times 5x1). Intermediate value
    if h2i > h2
        TH = Ti; % new right boundary
    else
        TL = Ti; % new left boundary
    end
end
T2 = (TH+TL)/2;
T2bis = T2;
for i=1:NSp
    hi2(i)    = HNasa(T2,SpS(i));
    si1(i)    = SNasa(T1,SpS(i));
    si2(i)    = SNasa(T2,SpS(i));
end
s1thermal = Yair*si1';
s2thermal = Yair*si2';
lnPr = (s2thermal-s1thermal)/Rg;                                            % ln(P2/P1) = (s2-s1)/Rg , see lecture (s2 are only the temperature integral)
Pr = exp(lnPr);
P2 = P1*Pr;
S1  = s1thermal - Rg*log(P1/Pref);                                          % Total entropy stage 1
S2  = s2thermal - Rg*log(P2/Pref);                                          % Total entropy stage 2
% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,1,2);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T1,T2);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P1/kPa,P2/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v1,v2);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h1/kJ,h2/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S1/kJ,S2/kJ);
%% Difference between two approaches: so close but not identical
fprintf('----------------------------------------------\n%8s| %9.4f %9.4f  [K]\n----------------------------------------------\n','T2-int vs T2-bis',T2int,T2bis);
%% Here starts your part (compressor,combustor,turbine and nozzle). ...
%% [2-3] Compressor                   
sPart = 'Compressor';

P3 = P2 * P3overP2;

% Thermal part of entropy at state 2
s2T = interp1(TR, sair_a, T2);

% Isentropic compression: S3 = S2
s3T = s2T + Rg*log(P3/P2);

% Find T3 from thermal entropy
T3 = interp1(sair_a, TR, s3T);

% Find h3 from T3
h3 = interp1(TR, hair_a, T3);

% Velocity change neglected in compressor
v3 = 0;

Qdot_c = 0;% Adiabatic compressor

% Compressor specific work(assuming v3=0)
wc = h3 - h2 + 0.5*(v3^2-v2^2);

mair = AF * mfurate;

% Compressor work input
Wdot_c = mair*wc - Qdot_c;


% Calculate properties directly at T3
for i = 1:NSp
    h3i(i) = HNasa(T3,SpS(i));
    s3i(i) = SNasa(T3,SpS(i));
end

h3check   = Yair*h3i';
s3thermal = Yair*s3i';

% Total entropy
S2 = s2T       - Rg*log(P2/Pref);
S3 = s3thermal - Rg*log(P3/Pref);

fprintf('\n');
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,2,3);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T2,T3);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P2/kPa,P3/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v2,v3);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h2/kJ,h3/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S2/kJ,S3/kJ);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f  [kJ/kg]\n','wc',wc/kJ);
fprintf('%8s| %9.2f  [kW]\n','Wdot_c',Wdot_c/kJ);
%% [3-4] Combustor
%----MASS----
sPart = 'Combustor';

mf    = mfurate;          % Fuel mass flow [kg/s]
mair  = AF * mf;          % Air mass flow [kg/s]
mdot4 = mair + mf;        % Product mass flow [kg/s]

% Xair = [H2 O2 CO2 H2O N2]

r = Xair(5)/Xair(2);     % N2/O2 molar ratio 


a = AF*Mi(1) / (Mi(2) + r*Mi(5));% convert the  AF into the molar coefficient a


% Inlet mole amounts
N3 = [1, a, 0, 0, r*a];%compostion at the inlet [H2 O2 CO2 H2O N2]

% Inlet mole fractions
X3 = N3 / sum(N3);%compostion into ratio

% Convert mole amounts to corresponding masses
M3 = N3 .* Mi;

% Inlet mass fractions
Y3 = M3 / sum(M3);

N4 = [0, a-0.5, 0, 1, r*a];%composition at the outlet [H2 O2 CO2 H2O N2]

% Product mole fractions
X4 = N4 / sum(N4);%composition out ratio

% Convert mole amounts to corresponding masses
M4 = N4 .* Mi;

% Product mass fractions
Y4 = M4 / sum(M4);

%----Thermo computations----
% Isobaric combustor
P4 = P3;


v4 = v3;% Neglect velocity change through combustor

Qdot = 0;% Adiabatic combustor
Wdot = 0;% No work

Tfuel = Tamb;% Fuel inlet temperature

hf = HNasa(Tfuel,SpS(1));% Specific enthalpy of H2 fuel



h4 = (Qdot - Wdot + mair*h3 + mf*hf) / mdot4;% Energy conservation at v3=v4=0 W=0 Q=0

h4_a = Y4 * hia';


T4 = interp1(h4_a, TR, h4);% Find temperature corresponding to h4


S4 = 0;

for i = 1:NSp

    if X4(i) > 0

        s4i = SNasa(T4,SpS(i));%temperature part of the species i
        Ri = Runiv/Mi(i);%specific gas const of the species i
        P4i = X4(i)*P4;%partial pressure of the species i

        S4 = S4 + Y4(i) *(s4i - Ri*log(P4i/Pref));

    end
end