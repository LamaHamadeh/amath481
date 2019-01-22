function rhs = vorticityBackslash(t, w_vec, rhsFun, A) 

psi_vec = A\w_vec;





rhs = rhsFun(psi_vec, w_vec);

end