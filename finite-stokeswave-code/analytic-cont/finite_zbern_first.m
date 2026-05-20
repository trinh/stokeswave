function zp = finite_zbern_first(f, z, wave)

mu = wave.mu;
B = wave.B;

[zminus, zpminus] = wave.getZValues(-f);

zp = mu./(2.*zpminus.*(mu*B + 1i*pi*(z + zminus)));