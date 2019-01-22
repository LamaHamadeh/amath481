function rhs=rhspde(t,u,k,dx,A)
rhs=(k/dx^2)*A*u;