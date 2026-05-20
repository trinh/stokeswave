function [D, L, UD] = region_indicator(P,r0)


Q = -log(r0)/(2*pi);

if Q==Inf
    
        for j=1:length(P)
        if isnan(P(j))==1
                D(j)=NaN;
        else
            if 0>imag(P(j))
                D(j)= 1;
            else
                D(j)=0;
            end
        end
        end
else

        for j=1:length(P)
            if isnan(P(j))==1
                    D(j)=NaN;
            else
                for n=0:51
                    if (-2*(n-1)*Q>imag(P(j)) && imag(P(j))>=-2*(n)*Q) 
                        D(j)= n;
                end
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


UD =[];
for j=1:length(L)
    UD(j)=find_updown(j+1, L);
end

L = [0, L];
UD = [1, UD];


%P is a vector indicating the order of line-crossing (P(j)=n indicates
%the jth line crossed is -2Qni


 
     
   
    
    

            