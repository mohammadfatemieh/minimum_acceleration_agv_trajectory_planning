function cost_mat = cost(fourier_order, derivative_order, keyframe_list)
  keyframe_cnt = size(keyframe_list, 1);
  time_idx = 3;

  cost_mat = [];
  for i = 1 : keyframe_cnt - 1
    L = 2 * (keyframe_list(i + 1, time_idx) - keyframe_list(i, time_idx));
    [A, w, p, type_ind] = fourier_sin(L, fourier_order);
    for i = 1 : derivative_order
      [A, w, p, type_ind] = fourier_der(A, w, type_ind);
    endfor
    cost_mat = blkdiag(cost_mat, ...
      fourier_square_int_val(A, w, type_ind, keyframe_list(i + 1, time_idx)), ...
      fourier_square_int_val(A, w, type_ind, keyframe_list(i, time_idx)))
  endfor
endfunction

function int_p = polyint(p, order)
  int_p = p;
  for i = 1 : min(length(p), order + 1)
    int_p(i) /= (order + 2 - i);
  endfor
endfunction

function [int_A, int_w, int_type_ind] = fourier_int(A, w, type_ind)
  int_A = A ./ w;
  int_w = w;
  int_type_ind = type_ind ./ 1i;
endfunction

function [val] = fourier_square_int_val(A, w, type_ind, x)
  n = size(A, 1);
  val = zeros(n, n);
  for i = 1 : n
    val(i, i) = - A(i)^2 * (sin(2 * w(i) * x) - 2 * w(i) * x) / (4 * w(i));
  endfor
endfunction
