function zp = finite_zbern_flipode(s, z, wave, Zminus, Zpminus)
    % Operates on arbitrary parameter s 
    % Zminus and Zpminus take in s parameter expecting coordinate change
    % gam(s)

mu = wave.mu;
B = wave.B;

zminus = Zminus(s);
zpminus = Zpminus(s);

zp = mu./(2.*zpminus.*(mu*B + 1i*pi*(z + zminus)));

