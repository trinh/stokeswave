function G = getint_gauss(w, solseries)
% GETINT_GAUSS computes beta integral using a specialized scheme

F = @(u) imag(zpfree(u, solseries)).*cot(pi*(u - w));

I = pvi_gauss(F, w, -0.5, 0.5);

I = reshape(I, size(w));
G = 1 - I;

end