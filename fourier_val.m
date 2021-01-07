function [val] = fourier_val(A, w, p, type_ind, x)
  if all(real(type_ind)) ~= 0
    val = A .* sin(w .* x .- p) .* type_ind;
  else
    val = A .* sin(w .* x .- p) .* type_ind ./ i;
  endif
endfunction
