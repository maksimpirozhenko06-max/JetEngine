clear all;
close all;
clc;
global Runiv Pref
Runiv = 8.314472;
Pref  = 1.01235e5;

abspath_to_generalfolder = ...
    "C:\Users\20251791\Downloads\JetEngineMatlab\JetEngineMatlab\General";

addpath(abspath_to_generalfolder);

TdataBase = fullfile(abspath_to_generalfolder, 'NasaThermalDatabase');
load(TdataBase);

whos

O2ind = myfind({Sp.Name}, {'O2'});
Oind = myfind({Sp.Name}, {'O'});
Nind = myfind({Sp.Name}, {'N'});
T = 300:10:4000;
T0 = 25+273.15;

hO2 = HNasa(T0, Sp(O2ind))
hO1 = HNasa(T0,Sp(Oind))
N = HNasa(T0,Sp(Nind))
%Exercise 2
N2ind = myfind({Sp.Name}, {'N2'});
hN2=HNasa(T,Sp(N2ind));
hO2 = HNasa(T, Sp(O2ind));
MO2=0.016*2
MN2=0.014*2
YN2 = (0.79*MN2) / (0.79*MN2 + 0.21*MO2);
YO2 = (0.21*MO2) / (0.79*MN2 + 0.21*MO2);

hO2N2 = hN2*YN2 + hO2*YO2;
plot(hO2N2,T);