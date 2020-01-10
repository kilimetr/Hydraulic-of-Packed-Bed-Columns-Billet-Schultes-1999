close all; clear all; clc;

% PHASES CHARACTERISTICS
%B           = 106.1825;      % liq velocity     [m3/m2/hod]
etaL        = 0.001021;      % liq dynam viscos [Pas]
rhoL        = 998.773;       % liq density      [kg/m3]
uV          = 1.5;           % gas velocity     [m/s]
etaV        = 1.8*10^(-5);   % gas dynam viscos [Pas]
rhoV        = 1.23595;       % gas density      [kg/m3]
g           = 9.81;          % grav accelerat   [m/s2]

% PACKING CHARACTERISTICS
epsilon     = 0.945;           % void fraction         [-]     RALU PAK YC-250
a           = 250;             % specific surface area [m2/m3]
CFltab      = 2.558;           % packings constant     [-]
CStab       = 3.178;           % packings constant     [-]
Cp0         = 0.191;           % packings constant     [-]
CL          = 1.334;           % packings constant     [-]
CV          = 0.385;           % packings constant     [-]
dp          = 6*(1-epsilon)/a; % particle diameter     [m]

% COLUMN CHARACTERISTICS
dkol = 1;               % column diameter [m]
Skol = pi * dkol^2 / 4; % column surface  [m2]

VFl = rhoV*uV*Skol*3600; % gas mass flow rate [kg/s]

% CALCULATIONS - Flooding Point
uLFl  = 1;
hLmin = epsilon/3;
hLmax = epsilon;
hLFl  = (hLmin + hLmax)/2;
psi   = 1;

pars = [uV g epsilon a rhoL rhoV VFl CFltab Skol etaL etaV];

result = fsolve(@(y) BS_fl(pars,y),[uLFl hLFl psi]);
disp(result); disp(BS_fl(pars,result));

uLFl = result(1);
hLFl = result(2);
psi  = result(3);

% CALCULATIONS - Loading Point
uLFl = 0.0295;
uLS = uLFl;

LS = rhoL*uLS*Skol*3600;
parss = [uLS g epsilon a rhoL rhoV LS CStab Skol etaL etaV];

uVS  = 5;
psiS = 5;

result2 = fsolve(@(y) BS_ld(parss,y),[uVS psiS]);
disp(result2); disp(BS_ld(parss,result2));

uVS  = result2(1);
psiS = result2(2);
uL   = uLFl;
hLS  = (12 * 1/g * etaL/rhoL * uL * a^2)^(1/3); % ok

% DRY PRESSURE DROP
uV   = 1.5;
ds   = dkol;
K    = 1 / (1 + 2/3 * 1/(1-epsilon) * dp/ds);
ReV  = uV*dp / ((1-epsilon)*etaV) * rhoV*K;
psi0 = Cp0 * (64/ReV + 1.8/(ReV^(0.08)));
FV   = uV * sqrt(rhoV); % gas load factor
dp0H = psi0 * a/(epsilon^3) * (FV^2)/2 * 1/K;

% WET PRESSURE DROP
uVFl = uV;
hL   = hLS + (hLFl - hLS)*((uV/uVFl)^(13)); % uVS<uV<uVFl
C1   = 13300 / (a^(3/2));
FrL  = uL^2 * a/g;
psiL = Cp0 * (64/ReV + 1.8/(ReV^0.08)) * ((epsilon-hL)/epsilon)^(1.5) * (hL/hLS)^(0.3) * ...
       exp(C1*sqrt(FrL));
dpH  = psiL * a/((epsilon-hL)^3) * (FV^2)/2 * 1/K;




% COEFFICIENTS - for uV <= uVS
uLi   = uL/hLS;      % effective liquid velocity [m/s]
dh    = 4*epsilon/a; % hydraulic diameter [m]
aPha  = 1.5 * (a*dh)^(-0.5) * (uL*dh*rhoL/etaL)^(-0.2) * ...
        ((uL^2)*rhoL*dh/sigmaL)^(0.75) * ((uL^2)/(g*dh))^(-0.45); % [-]
kLaPh = CL*(12^(1/6)) * (uLi^0.5) * ((DL/dh)^0.5) * a *aPha; % volumetric mass transfer coeff [1/s]
kVaPh = CV/((epsilon-hLS)^0.5) * (a^(3/2))/(dh^0.5) * DV * ((uV*rhoV/(a*etaV))^(3/4)) * ...
        (etaV/(rhoV*DV))^(1/3) * aPha; % volumetric mass transfer coeff [1/s]
kxaPh = kLaPh * cL; % mass transfer coeff [mol/m3/s] recalculation
kyaPh = kVaPh * cV; % mass transfer coeff [mol/m3/s] recalculation
Kya   = 1 / (smernice/kxaPh + 1/kyaPh); % 
HETP  = log(smernice)/(smernice-1) * cV*uV/Kya; % 


