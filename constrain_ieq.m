function [mat, vec] = constrain_ieq(fourier_order, keyframe_list, vel_bound)
  fourier_size = 2 * fourier_order + 1;
  keyframe_cnt = size(keyframe_list, 1);
  time_idx = 3;
  mat = [];
  for i = 1 : keyframe_cnt - 1
    blk = [ones(1, fourier_size); -ones(1, fourier_size)];
    blk(:, 1) = 0;
    mat = blkdiag(mat, blk);
  endfor
  mat = blkdiag(mat, mat);
  vec = ones(4 * (keyframe_cnt - 1), 1);
endfunction
