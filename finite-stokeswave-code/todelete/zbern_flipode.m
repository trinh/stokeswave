function zp = zbern_flipode(t, z, solseries, Zminus, Zpminus)

mu = solseries.mu;
B = solseries.B;

zminus = Zminus(t);
zpminus = Zpminus(t);
zp = mu./(2.*zpminus.*(mu*B + 1i*pi*(z + zminus)));
