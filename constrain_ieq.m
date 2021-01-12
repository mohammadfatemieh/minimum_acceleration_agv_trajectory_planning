function [mat, vec] = constrain_ieq(fourier_order, keyframe_list, vel_bound)
  fourier_size = 2 * fourier_order + 1;
  keyframe_cnt = size(keyframe_list, 1);
  time_idx = 3;
  mat = [];
  for i = 1 : keyframe_cnt - 1
    pos_f = fourier(fourier_order, [keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx)]);
    vel_A = pos_f.derivative.amplitude;
    % blk = [diag(vel_A.^2); diag(-(vel_A.^2))]
    blk = [vel_A.^2; -vel_A.^2];
    % blk = [diag(ones(1, fourier_size)); diag(-ones(1, fourier_size))];
    blk(:, 1) = 0;
    mat = blkdiag(mat, blk);
  endfor
  % mat = blkdiag(mat, mat);
  mat = [mat, mat];
  %vec = ones(4 * fourier_size * (keyframe_cnt - 1), 1);
  %vec = ones(2 * fourier_size * (keyframe_cnt - 1), 1);
  vec = 0.71 * ones(2 * (keyframe_cnt - 1), 1);
  % vec = ones(4 * (keyframe_cnt - 1), 1);
endfunction
