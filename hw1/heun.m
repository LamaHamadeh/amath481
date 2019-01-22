function y = heun(F, y0, t_int, dt)

t = t_int(1):dt:t_int(2);

y = zeros(length(t) - 1, 1);
y(1) = y0;

for n = 1:length(t)-1
   
    intermediate = F(t(n) + dt, y(n) + dt*F(t(n), y(n)));
    y(n + 1) = y(n) + (dt/2)*(F(t(n), y(n)) + intermediate);
    
end
    
end