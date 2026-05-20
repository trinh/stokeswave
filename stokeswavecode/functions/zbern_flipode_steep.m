function zp = zbern_flipode_steep(t, z, steepsolseries, Zminus, Zpminus)
% ZBERN_FLIPODE calculates the dz/df values used in the rectangular
% flipping scheme
%
% Zminus = z(-f) 
% Zpminus = z'(-f)
% Solseries - truncation method data


mu = steepsolseries.mu;

L=10;
[~,zp] = zfree_steep(0, steepsolseries, L);

B = 0.5*(zp).^(-2);

zminus = Zminus(t);
zpminus = Zpminus(t);
zp = mu./(2.*zpminus.*(mu*B + 1i*pi*(z + zminus)));
