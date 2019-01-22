function y = forwardEuler(F, y0, t_int, dt)

t = t_int(1):dt:t_int(2);

y = zeros(length(t) - 1, 1);
y(1) = y0;

for n = 1:length(t)-1
   
    y(n + 1) = y(n) + dt * F(t(n), y(n));
    
end
    
end