function [A, w, p, type_ind] = fourier_sin(l, r, f_order)
  L = r - l;
  A = ones(f_order, 1) .* sqrt(2 / L);
  w = zeros(f_order, 1);
  p = zeros(f_order, 1);
  for i = 1 : f_order
    w(i) = pi * i / L;
    p(i) = -l * pi * i / L;
  endfor
  type_ind = ones(f_order, 1) * (1 + 0i);
endfunction
