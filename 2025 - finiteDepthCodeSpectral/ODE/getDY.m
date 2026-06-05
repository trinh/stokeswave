function [outP] = getDY(Y,F,C,Z,Ifunc)

    %Iz = Ifunc(Z);
    %N = length(Yreal);
    Iz = Ifunc(Z);
    if sign(imag(Z))==0
        a = 1; %Ensuring that we take the intial residue correctly 
    else
        a = sign(imag(Z));
    end
    outP = -a*1i/2.*((F^2/2)./((F^2/2+C-Y).*(1-Iz))-(1-Iz));
    if isnan(real(outP))
        %keyboard
    end
end
