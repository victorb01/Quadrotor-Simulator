%Function that outputs omegaDot
function omegaDot = omegaDot(t, omega, e_a, S)
    omegaDot = (S.quadParams.cm.*e_a - omega)./S.quadParams.taum;
end