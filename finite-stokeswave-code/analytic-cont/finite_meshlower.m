function [Z,Zp, F] = finite_meshlower(wave, lastZ, lastZp, F)

phi = real(F(1,:));
psi = imag(F(:,1));
    
ep = wave.ep;
r0 = wave.r0;
name = ['../../finite-stokeswave-data/recdata_ep', num2str(ep), '_r0', num2str(r0), '.mat'];
tmp = load(name); %load in initial curve data pts from rectangular scheme
ff0 = tmp.recdat(2).ff{4};
zz0 = tmp.recdat(2).zz{4};

for j = 1:length(phi)
    current_phi = phi(j);
    
    startz0 = interp1(ff0, zz0, current_phi); %construct initial curve
    
    lastZfn = @(s) interp1(psi, lastZ(:,end-j+1), imag(s)); 
    lastZpfn = @(s)interp1(psi, lastZp(:,end-j+1), imag(s));
    fwd = @(psi, z) 1i*finite_zbern_flipode(current_phi + 1i*psi, z, wave, lastZfn, lastZpfn);
    
    %solve along this line
    [psipsi, zz] = ode45(fwd, psi, startz0);
    
    %fill column j of Z matrix with these solved values
    Z(:, j) = zz; 
    Zp(:,j) = finite_zbern_flipode(current_phi + 1i*psipsi, zz, wave, lastZfn, lastZpfn);
end
