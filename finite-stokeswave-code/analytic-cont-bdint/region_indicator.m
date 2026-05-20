function [D,L] = region_indicator(P,r0)

Q = -log(r0)/(2*pi);

for j=1:length(P)
    if isnan(P(j))==1
            D(j)=NaN;
    else
    for n=1:51
        if (-2*(n-1)*Q>imag(P(j)) && imag(P(j))>=-2*(n)*Q) 
            D(j)= n;
        end
    end
    end
end

numpass=0;
P=[];
for j=1:length(D)
    if (isnan(D(j))==1)
        numpass=numpass+1;
        L(numpass) = min(D(j-1),D(j+1));
    end
end

%P is a vector indicating the order of line-crossing (P(j)=n indicates
%the jth line crossed is -2Qni


 
     
   
    
    

            