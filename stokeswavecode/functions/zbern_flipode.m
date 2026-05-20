function [zp, P] = zbern_flipode(t, z, solseries, Zminus, Zpminus)
% ZBERN_FLIPODE calculates the dz/df values used in the rectangular
% flipping scheme
%
% Zminus = z(-f) 
% Zpminus = z'(-f)
% Solseries - truncation method data

mu = solseries.mu;
B = solseries.B;

zminus = Zminus(t);
zpminus = Zpminus(t);
P = mu*B + 1i*pi*(z + zminus);
zp = mu./(2.*zpminus.*P);
