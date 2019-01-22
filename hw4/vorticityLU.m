function rhs = vorticityLU(t, w_vec, rhsFun, L, U, P) 

y = L\(P*w_vec);
psi_vec = U\y;




rhs = rhsFun(psi_vec, w_vec);

end