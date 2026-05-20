function zp = finite_zbern_flipode(s, z, wave, Zminus, Zpminus)

mu = wave.mu;
B = wave.B;

zminus = Zminus(s);
zpminus = Zpminus(s);

zp = mu./(2.*zpminus.*(mu*B + 1i*pi*(z + zminus)));