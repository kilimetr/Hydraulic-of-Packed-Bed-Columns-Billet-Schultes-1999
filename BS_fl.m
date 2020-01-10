function [res] = BS_fl(pars, yvec)
% FLOADING POINT

uVFl    = pars(1);
g       = pars(2);
epsilon = pars(3);
a       = pars(4);
rhoL    = pars(5);
rhoV    = pars(6);
VFl     = pars(7);
CFltab  = pars(8);
Skol    = pars(9);
etaL    = pars(10);
etaV    = pars(11);


uLFl = yvec(1); %disp(uLFl);
hLFl = yvec(2); %disp(hLFl);
psi  = yvec(3); %disp(psi); pause(1);


LFl = uLFl*rhoL*Skol*3600;

if LFl/VFl*sqrt(rhoV/rhoL) < 0.4
    nFl = -0.194;
    CFl = CFltab;

elseif LFl/VFl*sqrt(rhoV/rhoL) > 0.4
    nFl = -0.708;
    CFl = 0.6244*CFltab*((etaL/etaV)^0.1028);
end


res(1) = uVFl - (sqrt(2) * (g/psi)^(0.5) * (epsilon-hLFl)^(3/2) / (sqrt(epsilon)) * ...
         (hLFl/a)^(0.5) * sqrt(rhoL/rhoV));
res(2) = psi - (g/(CFl^2) * (LFl/VFl * sqrt(rhoV/rhoL) * (etaL/etaV)^(0.2))^(-2*nFl));
res(3) = hLFl^3 * (3*hLFl-epsilon) - (6/g * a^2 * epsilon * (etaL/rhoL) * ...
         (LFl/VFl) * (rhoV/rhoL) * uVFl);
% res(4) = uLFl - (rhoV/rhoL * LFl/VFl * uVFl);
res = res';

condition = LFl/VFl*sqrt(rhoV/rhoL); %disp(condition); pause(1);

end

