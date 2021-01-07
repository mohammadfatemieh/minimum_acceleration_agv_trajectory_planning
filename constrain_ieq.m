function [mat, vec] = constrain_ieq(poly_order, keyframe_list, vel_bound)
  keyframe_cnt = size(keyframe_list, 1);
  time_idx = 3;
  polynomial = ones(1, poly_order);
  vel_poly = polyder(polynomial);
  mat = cell(1, 2);
  vec = ones(12 * (keyframe_cnt - 1), 1) * vel_bound(2);
  for i = 1 : keyframe_cnt - 1
    t_sample = linspace(keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx), 5);
    sub_mat = zeros(4, poly_order, 2, length(t_sample) - 2);
    for t = 2 : length(t_sample) - 1
      for j = 1 : 2
        sub_mat(j, 1 : poly_order - 1, 1, t - 1) =  poly_sub_t(vel_poly, t_sample(t));
      endfor
      for j = 3 : 4
        sub_mat(j, 1 : poly_order - 1, 1, t - 1) = -poly_sub_t(vel_poly, t_sample(t));
      endfor
      for j = 1 : 2 : 3
        sub_mat(j, 1 : poly_order - 1, 2, t - 1) =  poly_sub_t(vel_poly, t_sample(t));
      endfor
      for j = 2 : 2 : 4
        sub_mat(j, 1 : poly_order - 1, 2, t - 1) = -poly_sub_t(vel_poly, t_sample(t));
      endfor
    endfor
    % s_mat(:, :, 1) = [sub_mat(:, :, 1, 1); sub_mat(:, :, 1, 2)];
    % s_mat(:, :, 2) = [sub_mat(:, :, 2, 1); sub_mat(:, :, 2, 2)];
    for j = 1 : 2
      s_mat = [];
      for t = 1 : length(t_sample) - 2
        s_mat = [s_mat; sub_mat(:, :, j, t);];
      endfor
      mat{j} = blkdiag(mat{j}, s_mat);
    endfor
  endfor
  mat = [mat{1}, mat{2}];
endfunction

function substed_poly = poly_sub_t(polynomial, t)
  poly_order = length(polynomial);
  substed_poly = polynomial;
  for i = 1 : poly_order
    substed_poly(i) *= t^(poly_order - i);
  endfor
endfunction
