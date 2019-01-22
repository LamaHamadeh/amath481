function [x, iter] = bisection(f, a, b, tol, maxIter)

    x = zeros(maxIter, 1);
    iter = 1;
    x(iter) = (a + b) / 2;
    
    while abs(f(x(iter))) > tol && iter <= maxIter
        
        if f(a) * f(x(iter)) > 0 % same sign
            a = x(iter);
        else 
            b = x(iter);
        end  
        
        iter = iter + 1;
        
        x_next = (a + b) / 2; 
        x(iter) = x_next;
             
    end
    
    x = x(1:iter);

end