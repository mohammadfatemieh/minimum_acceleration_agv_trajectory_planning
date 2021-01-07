function cost_mat = cost(poly_order, derivative_order, keyframe_list)
  keyframe_cnt = size(keyframe_list, 1);
  time_idx = 3;

  polynomial = ones(1, poly_order);
  for i = 1 : derivative_order
    polynomial = polyder(polynomial);
  end

  cost_blk = zeros(poly_order);
  cost_blk(1 : poly_order - derivative_order, 1 : poly_order - derivative_order) = polynomial' * polynomial;
  for i = 1 : poly_order - derivative_order
    cost_blk(i, :) = polyint(cost_blk(i, :), (poly_order - derivative_order) * 2 - i - 1);
  end

  cost_mat = cost_subs_t(poly_order, derivative_order, cost_blk, keyframe_list(1, time_idx), keyframe_list(2, time_idx));
  for i = 2 : keyframe_cnt - 1
    cost_mat = blkdiag(cost_mat, cost_subs_t(poly_order, derivative_order, cost_blk, keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx)));
  endfor
  cost_mat = blkdiag(cost_mat, cost_mat);
endfunction

function int_p = polyint(p, order)
  int_p = p;
  for i = 1 : min(length(p), order + 1)
    int_p(i) /= (order + 2 - i);
  endfor
endfunction

function cost_blk = cost_subs_t(poly_order, derivative_order, temp_cost_blk, cur_t, nxt_t)
  cost_blk = temp_cost_blk;
  for i = 1 : (poly_order - derivative_order)
    for j = 1 : (poly_order - derivative_order)
      cost_blk(i, j) *= ...
        abs(nxt_t^((poly_order - derivative_order) * 2 + 1 - i - j) ...
          - cur_t^((poly_order - derivative_order) * 2 + 1 - i - j));
    endfor
  endfor
endfunction
