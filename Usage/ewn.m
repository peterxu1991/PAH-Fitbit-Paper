function [Y] = ewn(X, alpha)
Y =  zeros(length(X),1);
for i= 1:length(X)
    if i ==1
        
      Y(i) = X(i);
      
    else
        
        Y(i) = alpha*X(i) + (1-alpha)*Y(i-1);
    end


end

