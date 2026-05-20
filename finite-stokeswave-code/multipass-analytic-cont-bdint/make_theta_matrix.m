function T = make_theta_matrix(L, UD);



T=zeros(length(L), length(L));


    
for l=2:1:length(L)
   
    
    
    n=L(l); %we have just crossed the line imf = -2nQi

    %first check whether we are passing through from above or below
    %updown = ±1 if we pass going up/down

    updown = UD(l); 
   
    if n>1 %Here we don't need to worry about whether theta has an UHP or LHP primary equation
        %NEED TO CHECK IF N=1
        for j=1:l-1
            if L(j) == n-1 %finding the functions argument f+2(n-1)Qi for which this line is 'virtual fluid floor'
                    T(j,l)=-updown; %argument f+2Qni


                
            elseif L(j)==n %finding the functions argument f+2Qi for....
                            %which this line is 'virtual free surface'
                    T(j,l)=updown*2; % new fn argument f+2Qni
            end
        end
        
          %%%%%%%%%%%%
    
    elseif n==1 %Here we don't need to worry about whether theta has an UHP or LHP primary equation
        %NEED TO CHECK IF N=1
        for j=1:l-1
            if L(j)==n-1 %finding the functions argument f for which this line is fluid floor
                    if j==1 
                        upperlower = 1;
                    else
                        upperlower = UD(j);
                    end
                
                   T(j,l)=upperlower*updown*1i; 

                
            elseif L(j)==n %finding the functions argument f+2Qi for....
                            %which this line is 'virtual free surface'
                    T(j,l)=updown*2; 
            end
        end     
  
          
          
          
    elseif n == 0
            if l == 2 %i.e this is the first cross of a line
                  T(1,l)=2; 
            elseif l>2

                    for j=1:l-1

                         if L(j) == 0 %picking out functions of argument f
                                    if j==1 %the value of a for this sheet
                                        upperlower = 1;
                                    else
                                        upperlower = UD(j);
                                    end  
                            T(j,l)=-upperlower*updown*2;
                         end
                    end

            end
     end
end
end
