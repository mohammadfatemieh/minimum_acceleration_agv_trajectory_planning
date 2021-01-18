function result = verify(fourier_order, coef_list, keyframe_list, dangerous_region)
  fourier_size = 2 * fourier_order + 1;
  keyframe_cnt = size(keyframe_list, 1);
  time_idx = 3;
  result = logical(1);
  for i = 1 : keyframe_cnt - 1
    time_sample = linspace(keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx), 10);
    pos_x_f = fourier(fourier_order, [keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx)]);
    pos_x_f = pos_x_f.scale(coef_list((i - 1) * fourier_size + 1 : i * fourier_size, 1)');
    pos_y_f = fourier(fourier_order, [keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx)]);
    pos_y_f = pos_y_f.scale(coef_list((i - 1 + keyframe_cnt - 1) * fourier_size + 1 : (i + keyframe_cnt - 1) * fourier_size, 1)');
    pos_x = zeros(1, length(time_sample));
    pos_y = zeros(1, length(time_sample));
    for k = 1 : length(time_sample)
      pos_x(k) = sum(pos_x_f.value(time_sample(k)));
      pos_y(k) = sum(pos_y_f.value(time_sample(k)));
      if k >= 2
        trajectory_seg = [pos_x(k - 1), pos_y(k - 1); pos_x(k), pos_y(k)];
        intersect_ret = logical(0);
        intersect_ret = intersect_ret || is_intersected(trajectory_seg, [dangerous_region(1, :); dangerous_region(2, :)]);
        intersect_ret = intersect_ret || is_intersected(trajectory_seg, [dangerous_region(1, :); dangerous_region(3, :)]);
        intersect_ret = intersect_ret || is_intersected(trajectory_seg, [dangerous_region(2, :); dangerous_region(4, :)]);
        intersect_ret = intersect_ret || is_intersected(trajectory_seg, [dangerous_region(3, :); dangerous_region(4, :)]);
        if intersect_ret
          result = logical(0);
          break;
        endif
      endif
    endfor
  endfor
endfunction
