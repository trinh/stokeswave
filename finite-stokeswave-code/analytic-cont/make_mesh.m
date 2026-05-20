

clear
close all


r0 = 0.2;
ep = 0.6;
wave = VBStokesWave(r0, ep);


%path boundaries
phi0 = 0.5;
psi0 = 1;
phisamples = 50;
psisamples = 50;
phi = linspace(-phi0, phi0, phisamples);
psi = linspace(0, -psi0, psisamples);
numbsheets = 3; %number of sheets we want
[PHI, PSI] = meshgrid(phi, psi);

F = PHI + 1i*PSI;
first = @(v)v(1);
for rind = 1:numbsheets
    Zvals{rind} = [];
    Zpvals{rind} = [];
    if rind == 1
        %use calculated vals on first LHP
       
        [Z, Zp] = finite_meshlower_first(wave, F);
        Zvals{rind}=Z;
        Zpvals{rind}=Zp;        
     
     else
        if first(factor(rind)) == 2
            %get Z and Zp vals on an UHP mesh
            [Z, Zp, F] = finite_meshupper(wave, Zvals{rind-1}, Zpvals{rind-1}, F);
            Zvals{rind}=Z;
            Zpvals{rind}=Zp; 
            upperF = F;
            
        else
            %get Z and Zp vals on a LHP mesh  
            
            %make sure range of finite_zpath_rectangle data at least as big as mesh region 
            %(otherwise initial curves will not be long enough)

            [Z, Zp, F] = finite_meshlower(wave, Zvals{rind-1}, Zpvals{rind-1}, F);
            Zvals{rind}=Z;
            Zpvals{rind}=Zp;
    
        end
    end
end
    
%plot the imaginary surface
figure(1) 
for j = 1:numbsheets
    hold on
    if first(factor(j)) == 2
      p = surf(real(F), -imag(F), imag(Zvals{j}));
      set(p, 'EdgeColor', 'none'); 
    else 
        p = surf(real(F), imag(F), imag(Zvals{j}));
set(p, 'EdgeColor', 'none'); 
    end
end
xlabel('phi'); ylabel('psi'); zlabel('imag(z)');
hold off

%plot the real surface
figure(2); figshift
for j = 1:numbsheets
    hold on
    if first(factor(j)) == 2
      p = surf(real(F), -imag(F), real(Zvals{j}));
      set(p, 'EdgeColor', 'none'); 
    else 
        p = surf(real(F), imag(F), real(Zvals{j}));
set(p, 'EdgeColor', 'none'); 
    end
end
xlabel('phi'); ylabel('psi'); zlabel('real(z)');
hold off



%}