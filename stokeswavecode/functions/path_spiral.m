function path_out = path_spiral(path_in, N, no_loops, width)

path_out = path_in;
lst = path_in(end);

t = linspace(0,no_loops*pi,N);
x = real(lst)+width*t.*cos(t);
y = imag(lst)+width*t.*sin(t);


path_out = [path_in, x+1i*y, lst];

end

