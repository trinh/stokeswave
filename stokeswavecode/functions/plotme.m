function [h1, h2] = plotme(F, Z, col, N)

h1 = surf(real(F), imag(F), Z, 'FaceColor', col, 'EdgeColor', 'none');
h2 = wireframe(real(F), imag(F), Z, N);
