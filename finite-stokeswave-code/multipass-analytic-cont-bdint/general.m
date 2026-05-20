
n=L(number_of_pass-1); %we have just crossed the line imf = -2nQi

%first check whether we are passing through from above or below
%updown = ±1 if we pass going up/down

if number_of_pass == 2 %if first line cross (imf=0)
    updown = -1
else
    for i=number_of_pass-1:-1:1 %check for last time a different critical line was crossed
        if L(i)<n
            updown=(-1)^(number_of_pass-1-i) 
            continue
        elseif L(i)>n
            updown=(-1)^(number_of_pass-1-i+1)
            continue
        else %only an option if n=0
            updown = (-1)^(number_of_pass-1)
        end
    end
end


if n>0
    for j=1:number_of_pass-1
        if L(j)==n-1 %finding the functions argument f+2(n-1)Qi for which this line is 'virtual fluid floor'
            theta(j,:)=theta(j,:)-updown/1i*...
                (tau(number_of_pass,:)-1i*theta(number_of_pass,:))/2); %argument f+2Qni


        elseif L(j)=n %finding the functions argument f+2Qi for....
                        %which this line is 'virtual free surface'
                theta(j,:)=theta(j,:)-updown/1i*2*theta(number_of_pass,:); % new fn argument f+2Qni
        end
    end
      %%%%%%%%%%%%
elseif n=0
        if number_of_pass == 2 %i.e this is the first cross of a line
             theta(j,:)=theta(j,:)-updown/1i*2*theta(number_of_pass,:); 
        else

            for i=1:number_of_pass-1
                
                 if L(i) = 0 %picking out functions of argument f
                    for k=i:-1:1
%                         if L(k)>n %ie if theta(j,:) was introduced in UHP
%                              theta(j,:)=theta(j,:)+((-1)^updown)/1i*2*theta(number_of_pass,:);
%                         elseif
                             for i=j:-1:1 %check for last time a different critical line was crossed
                        if L(i)<0
                            upperlower=(-1)^(number_of_pass-1-i) 
                            continue
                        elseif L(i)>0
                            upperlower=(-1)^(number_of_pass-1-i+1)
                            continue
                        end
                             end
                    theta(j,:)=theta(j,:)+upperlower*((-1)^updown)/1i*2*theta(number_of_pass,:);
                    end
             
                     %GO BACK OVER THESE LAST FEW LINES
                     
                     
                     
                     
                     
                     
%                                      elseif L(j-1)<n %ie if theta(j,:) was introduced in LHP
%                                          theta(j,:)=theta(j,:)-((-1)^updown)/1i*2*theta(number_of_pass,:);
% 
%                                      elseif L(j-1)=n
%                                          %find first L(i), i<j st L(j-1)=/=n and repeat check.
%                                            for i=1:j-1
%                                                if L(i)=0
%                                                    continue
%                                                else %we know L(i)>0, j-i tells us whether theta(j,:) was UHP or LHP
%                                                    theta(j,:)=theta(j,:)+(-1)^(j-i)*...
%                                                        ((-1)^updown)/1i*2*theta(number_of_pass,:); 
%                                                end
%                                            end
%                                      end
%              end
%              
%              %Check whether theta(j,:) orginated in UHP or LHP
%              
%        
%     end
% end

       %         elseif L(j)=n+1 
%                          %find first L(i), i<j st L(j-1)=/=n and repeat check.
%                            for i=1:j-1
%                                if L(i)=0
%                                    continue
%                                else %we know L(i)>0, j-i tells us whether theta(j,:) was UHP or LHP
%                                    theta(j,:)=theta(j,:)+(-1)^(j-i)*...
%                                        ((-1)^updown)/1i*2*theta(number_of_pass,:); 
%                                end
% 
%                            end          
             
    