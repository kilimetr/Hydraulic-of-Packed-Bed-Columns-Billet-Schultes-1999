function [res] = BS_ld(pars, yvec)
% LOADING POINT

uLS     = pars(1);
g       = pars(2);
epsilon = pars(3);
a       = pars(4);
rhoL    = pars(5);
rhoV    = pars(6);
LS      = pars(7);
CStab   = pars(8);
Skol    = pars(9);
etaL    = pars(10);
etaV    = pars(11);

uVS  = yvec(1);
psiS = yvec(2);

VS = uVS*rhoV*Skol*3600;

if LS/VS*sqrt(rhoV/rhoL) < 0.4
    nS = -0.326;
    CS = CStab;
    
elseif LS/VS*sqrt(rhoV/rhoL) > 0.4
    nS = -0.723;
    CS = 0.695*CStab*(etaL/etaV)^(0.1588);
end
    

res(1) = uVS - (sqrt(g/psiS) * (epsilon/(a^(1/6)) - a^(1/2) * ...
         (12 * 1/g * etaL/rhoL * uLS)^(1/3)) * (12 * 1/g * etaL/rhoL * uLS)^(1/6) * ...
         sqrt(rhoL/rhoV));
res(2) = psiS - (g/(CS^2) * (LS/VS * sqrt(rhoV/rhoL) * (etaL/etaV)^(0.4))^(-2*nS));

res = res';

end

