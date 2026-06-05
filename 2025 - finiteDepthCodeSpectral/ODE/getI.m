function out = getI(Yreal,Z,N)
%Function to calculat I via spectral method
    
    [ny, nx] = size(Z); 
    out = zeros(ny,nx);
    L=0.5;
    dxi=2*L/N;
    for county = 1:ny
        ep = imag(Z(county,1));
        
        
        yR = Yreal;
        
        Mn = [0:N/2-1 0 -N/2+1:-1];
        n=[0:N/2-1 0 -N/2+1:-1];
        Mk=1i*n*pi/L;
        h = mean(yR);
        p = floor(ep/2*h); %p form notes
        
        %I am getting a bug when h is large we get weird errors with the
        %compution of exp(...). 
        temp = -2*pi*Mn*(ep-h*(1+2*p));
        temp = max(min(temp,400),-400);
        MhC = -1i.*exp(temp)./sinh(-2*pi*Mn*h);
        %MhC = -1i*sign(Mn).*(exp(abs(temp)))./(sinh(2*pi*h.*abs(Mn)));
        

        %K=0 case
        if sin(pi*ep/h)>=0
            MhC(1) = 1i*(-ep/h+1); 
            MhC(N/2+1) = 1i*(-ep/h+1);
        else
            MhC(1) = 1i*(-ep/h-1); 
            MhC(N/2+1) = 1i*(-ep/h-1);
        end
      
        HCyR1 = ifft(MhC.*Mk.*fft(yR));
        %keyboard
        out(county,:) = interp1(linspace(-0.5, 0.5-dxi, N),HCyR1,real(Z(county,:)));
   end
end

