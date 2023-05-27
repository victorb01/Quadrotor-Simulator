%Function that returns aerodynamic drag
function d_a = aeroD(z_I, v_I, P)
    d_a = 0.5*P.quadParams.Cd*P.quadParams.Ad*P.constants.rho*v_I*dot(z_I,v_I);
end