function [der_A, der_w, der_p, der_type_ind] = fourier_der(A, w, p, type_ind)
  der_A = A .* w;
  der_w = w;
  der_p = p;
  der_type_ind = type_ind .* 1i;
endfunction
