% Using rectangular data

%close all
%clear

% =====================================
%         Load in the solution
F = load('../../stokeswavedata/recdat_steep_rec2_ep0.9.mat');
recdat = F.recdat;

% F = load('../../stokeswavedata/series_N100_ep0.3.mat');
% solseries = F.solseries;
% =====================================


rind = 2; % 2 or greater
integral = 0;
figure(1);
hold on

% singularity location
for k=1:8
    if mod(k,4) == 1
        t = imag(recdat(rind).ff{k})*((-1)^(rind));
        dgam = ((-1)^(rind))*1i;
    elseif mod(k,4) == 2
        t = real(recdat(rind).ff{k})*((-1)^(rind));
        dgam = ((-1)^(rind));
    elseif mod(k,4) == 3
        t = -imag(recdat(rind).ff{k})*((-1)^(rind));
        dgam = -1i*((-1)^(rind));
    else 
        t = -real(recdat(rind).ff{k})*((-1)^(rind));
        dgam = -1*((-1)^(rind));
    end
    
   % gg = recdat(rind).gg{k} - 2i*solseries.B*solseries.mu/(2*pi);
   % ggp = recdat(rind).ggp{k};
    ff = recdat(rind).ff{k};
    zz = recdat(rind).zz{k};
    zzp = recdat(rind).zzp{k};
    zzpp = central_diff(zzp, t);
    
    integrand = (ff.*(zzpp./zzp));
    integral = integral + trapz(t, integrand);
   
    
end
loc = (1/(2*pi*1i))*integral;
disp(['singularity location, f = ', num2str(loc, '%15.12f')]);


%{

 % And the coefficient
 
 for k=1:8
    if mod(k,4) == 1
        t = imag(recdat(rind).ff{k})*((-1)^(rind));
        dgam = ((-1)^(rind))*1i;
    elseif mod(k,4) == 2
        t = real(recdat(rind).ff{k})*((-1)^(rind));
        dgam = ((-1)^(rind));
    elseif mod(k,4) == 3
        t = -imag(recdat(rind).ff{k})*((-1)^(rind));
        dgam = -1i*((-1)^(rind));
    else 
        t = -real(recdat(rind).ff{k})*((-1)^(rind));
        dgam = -1*((-1)^(rind));
    end
    
    ff = recdat(rind).ff{k};
    zz = recdat(rind).zz{k};
    zzp = recdat(rind).zzp{k};
    zzpp = central_diff(zzp, t);
    gg = recdat(rind).gg{k} - 2i*solseries.B*solseries.mu/(2*pi);
    ggp = recdat(rind).ggp{k};
    
    
    integrand = dgam.*gg.*(ff-loc).^(-3/2);
    integral = integral + trapz(t, integrand);
    
    figure(4);
    hold on
        plot3(real(ff), imag(ff), imag(zzp));
   
end

coeff = (1/(4*pi*1i))*integral

%}

