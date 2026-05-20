function [phi,psi,FP]=SeriesData(wave, pathdat)
            % Plot the solution
            wave=wave;
            psipath=imag(pathdat.f);
            phi = linspace(-0.5,0.5,40);
            psi = linspace(min(psipath), 0, 40);
            
            [PHI, PSI] = meshgrid(phi, psi);
            F = PHI + 1i*PSI;
            
            [Z, ZP] = getZValues(wave, F);
            
            FP= 1./ZP;
            
            meshgrid(phi,psi,FP);
            
            
            