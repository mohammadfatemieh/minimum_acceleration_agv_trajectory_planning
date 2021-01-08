function solution = solve(poly_order, der_order, keyframe_list, end_point_cond, vel_bound)
  keyframe_cnt = size(keyframe_list, 1);
  time_idx = 3;
  tolr = 2e-4;
  step_len = 0.1;
  der_len = 0.001;
  vel_cond = end_point_cond(1 : 2, 4 : 5);
  acc_cond = end_point_cond(1 : 2, 6 : 7);

  H = cost(poly_order, der_order, keyframe_list);
  f = zeros(2 * (keyframe_cnt - 1) * (2 * poly_order + 1), 1);

  if keyframe_cnt == 2
    [A, b] = constrain_eq(poly_order, keyframe_list, vel_cond, acc_cond);
    [solution, _] = quadprog(H, f, [], [], A, b);
    plot_traj(poly_order, poly_coef, keyframe_list);
    return;
  endif

  prev_value = 127;
  cur_value = -128;

  dir_vec = -ones(keyframe_cnt - 1, keyframe_cnt) / (keyframe_cnt - 2);
  dir_vec(:, 1) = 0;
  for i = 1 : keyframe_cnt - 1
    dir_vec(i, i + 1) = 1;
  endfor
  dir_vec_len = sqrt(1 + 1 / (keyframe_cnt - 2));
  dir_vec = dir_vec' / dir_vec_len;

  while abs((prev_value - cur_value) / prev_value) > tolr
    step_len = (prev_value - cur_value) / prev_value * 0.1;
    prev_value = cur_value;
    vel_cond
    acc_cond
    [A, b] = constrain_eq(poly_order, keyframe_list, vel_cond, acc_cond);
    tic;
    [poly_coef, cur_value] = quadprog(H, f, [], [], A, b);
    toc;
    % DEBUG
    cur_value
    % DEBUG
    % plot_traj(poly_order, poly_coef, keyframe_list);
    grad = zeros(keyframe_cnt, 1);
    for i = 1 : keyframe_cnt - 1
      nkl = new_keyframe_list(keyframe_list, dir_vec(:, i), der_len, i);
      [A, b] = constrain_eq(poly_order, nkl, vel_cond, acc_cond);
      tic;
      [poly_coef, nxt_value] = quadprog(H, f, [], [], A, b);
      toc;
      grad += dir_vec(:, i) * (cur_value - nxt_value) / der_len;
    endfor
    keyframe_list(:, time_idx) += prefix_sum(step_len * grad);
  endwhile
  plot_traj(poly_order, poly_coef, keyframe_list);
  solution = poly_coef;
endfunction

function sum_array = prefix_sum(array)
  sum_array = array;
  for i = 2 : length(array)
    sum_array(i) += sum_array(i - 1);
  endfor
endfunction

function new = new_keyframe_list(old, dir_vec, der_len, i)
  keyframe_cnt = size(old, 1);
  new = old;
  new(:, 3) += der_len * prefix_sum(dir_vec);
endfunction
