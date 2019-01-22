function [x, iter] = newtonRaphson(f, df, x1, tol, maxIter)

    x = zeros(maxIter, 1);
    iter = 1;
    x(iter) = x1;
    
    while abs(f(x(iter))) > tol && iter <= maxIter
        
        x_next = x(iter) - f(x(iter))/df(x(iter));
        
        iter = iter + 1;
        x(iter) = x_next;
        
    end
    
    x = x(1:iter);
end