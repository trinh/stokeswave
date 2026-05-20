function A = split_path(path, r0)

A = path;
pm = 1i*1e-5;
Q = -log(r0)/(2*pi);


for n=0:50 %large enough for any path 50 times depth
    b=-2*Q*n;
    for j=1:100000 %ensures loop is longer than as yet unknown size of A
        if j<= length(A)-1
            if (imag(A(j))< b) && (imag(A(j+1))> b)
                ratio=(b-imag(A(j)))/(imag(A(j+1))-imag(A(j)));
                
                a = real(A(j))+ratio*(real(A(j+1))-real(A(j))); 
                
                   BEFORE = A(1:j);
                   AFTER = A(j+1:end);
                   A=[BEFORE a+b*1i-pm NaN*(1+1i) a+b*1i+pm AFTER];
                  
a = (b-imag(A(j)))/(imag(A(j+1))-imag(A(j)))*(real(A(j+1))-real(A(j)));

            elseif (imag(A(j))> b) && (imag(A(j+1))< b)
                    ratio=(b-imag(A(j)))/(imag(A(j+1))-imag(A(j)));
                    a = real(A(j))+ratio*(real(A(j+1))-real(A(j))); 
                
                   BEFORE = A(1:j);
                   AFTER = A(j+1:end);
                   A=[BEFORE a+b*1i+pm NaN*(1+1i) a+b*1i-pm AFTER];
            else
                continue
            end
        end
    end
end


            
           
        
        
    
  