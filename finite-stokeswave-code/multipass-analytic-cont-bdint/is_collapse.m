function [collapse,x0] = is_collapse(x0, collapsetol, Q)

C=zeros(1,50);
EXOS=zeros(1,50);
for n=1:50
    if abs(imag(x0)+2*n*Q) < collapsetol
        EXOS(n) = real(x0)-2*n*Q*1i;
        C(n) = 1;

    else
        C(n) = 0;
        EXOS(n)=NaN;

    end
end


collapse = max(C);

for j=1:50
    if isnan(EXOS(j))==0
        x0= EXOS(j);
    end
end