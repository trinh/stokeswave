function G = qt_getint(myz, solqtbdint)

phi = solqtbdint.phi;
theta = solqtbdint.theta;

N = solqtbdint.N;
dphi = phi(2) - phi(1);

getz = myz(:);
I = 0*getz;

for s = 1:length(getz)
    g = 2*theta.*cot(pi*(phi - getz(s)));
    g(1) = g(1)/2;
    g(N) = g(N)/2;
    I(s) = dphi/2*sum(g);
end

I = reshape(I, size(myz));
G = I; % This is essentially logq

end