function updown = find_updown(sheetnumber, L)


l=sheetnumber;
n=L(l-1);

if l ==1
    warning('Using find_updown without any line crossing');    
end

if l == 2 %if first line cross (imf=0)
        updown = -1;
        return
end

for i=l-1:-1:1 %check for last time a different critical line was crossed
            if L(i)<n
                updown=(-1)^(l-1-i);
                return
            end
            if L(i)>n
                updown=(-1)^(l-1-i+1);
                return
            end
            if i==1;%only an option if n=0
                updown = (-1)^(l-1);
                return
            end
 
        end
end
