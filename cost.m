function cost_mat = cost(fourier_order, derivative_order, keyframe_list)
  keyframe_cnt = size(keyframe_list, 1);
  time_idx = 3;

  cost_mat = [];
  for i = 1 : keyframe_cnt - 1
    f = fourier(fourier_order, [keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx)]);
    for _ = 1 : derivative_order
      f = f.derivative;
    endfor
    % Square the Fourier series and then integrate
    cost_mat = blkdiag(cost_mat, ...
      diag(f.square_int_value( ...
        keyframe_list(i, time_idx), ...
        keyframe_list(i + 1, time_idx) ...
      )) ...
    );
  endfor
  cost_mat = blkdiag(cost_mat, cost_mat);
endfunction
