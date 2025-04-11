function u_crt = complex_correction(u,u_ref)
u_crt = u/ (sum(u.*conj(u_ref),"all")/sum(abs(u_ref).^2,"all"));
end
