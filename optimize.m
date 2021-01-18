function [solution, keyframe_list] = optimize(fourier_order, der_order, keyframe_list, end_point_cond, vel_bound)
  % optimize  Optimize the time distribution for each segment.
  %
  %   [solution, keyframe_list] = optimize(fourier_order, der_order, keyframe_list, vel_bound)
  %     Return the optimal time distribution that minimize the target function
  %     (see bottom of optimize.m)
  keyframe_cnt = size(keyframe_list, 1);
  time_idx = 3;
  tolr = 1e-3;
  step_len = 0.0000001;
  der_len = 0.000001;
  vel_cond = end_point_cond(1 : 2, 4 : 5);
  acc_cond = end_point_cond(1 : 2, 6 : 7);
  m = 20;
  s = zeros(keyframe_cnt, m);
  g = zeros(keyframe_cnt, m);
  y = zeros(keyframe_cnt, m);
  rho = zeros(1, m);

  f = zeros(2 * (keyframe_cnt - 1) * (2 * fourier_order + 1), 1);

  prev_value = 127;
  cur_value = -128;

  [solution, fval] = solve(fourier_order, der_order, keyframe_list, end_point_cond, vel_bound);
  cur_value = target(fourier_order, solution, fval, keyframe_list, vel_bound);

  if keyframe_cnt == 2
    return;
  endif

  g(:, 1) = gradient(fourier_order, der_order, der_len, keyframe_list, end_point_cond, vel_bound, cur_value);
  s(:, 1) = step_len * g(:, 1);
  keyframe_list(:, time_idx) += prefix_sum(s(:, 1));
  prev_value = cur_value;
  [solution, cur_value] = solve(fourier_order, der_order, keyframe_list, end_point_cond, vel_bound);
  [min_vel, max_vel] = velocity_range(fourier_order, solution, keyframe_list);
  cur_value = (min_vel - vel_bound(1))^2 + (max_vel - vel_bound(2))^2;
  g(:, 2) = gradient(fourier_order, der_order, der_len, keyframe_list, end_point_cond, vel_bound, cur_value);
  y(:, 1) = g(:, 2) - g(:, 1);
  rho(1) = 1 / (y(:, 1)' * s(:, 1));

  % quasi-Newton method, BFGS algorithm
  for k = 2 : m
    alpha = zeros(1, m);
    beta = zeros(1, m);
    q = g(:, k);
    for i = min(m, k - 1) : -1 : 1
      alpha(i) = rho(i) * s(:, i)' * q;
      q = q - alpha(i) * y(:, i);
    endfor
    H0 = (s(:, k - 1)' * y(:, k - 1)) / (y(:, k - 1)' * y(:, k - 1)) * diag(ones(1, keyframe_cnt));
    z = H0 * q;
    for i = 1 : max(m, k - 1)
      beta(i) = rho(i) * y(:, i)' * z;
      z = z + s(:, i) * (alpha(i) - beta(i));
    endfor
    z = -z;
    s(:, k) = z;
    keyframe_list(:, time_idx) += prefix_sum(s(:, k));
    prev_value = cur_value;
    [solution, fval] = solve(fourier_order, der_order, keyframe_list, end_point_cond, vel_bound);
    cur_value = target(fourier_order, solution, fval, keyframe_list, vel_bound);
    cur_value
    if (prev_value - cur_value) / prev_value <= tolr
      break;
    endif
    g(:, k + 1) = gradient(fourier_order, der_order, der_len, keyframe_list, end_point_cond, vel_bound, cur_value);
    y(:, k) = g(:, k + 1) - g(:, k);
    rho(k) = 1 / (y(:, k)' * s(:, k));
  endfor

  plot_traj(fourier_order, solution, keyframe_list);
endfunction

function grad = gradient(fourier_order, der_order, der_len, keyframe_list, end_point_cond, vel_bound, cur_value)
  keyframe_cnt = size(keyframe_list, 1);
  time_idx = 3;

  dir_vec = -ones(keyframe_cnt - 1, keyframe_cnt) / (keyframe_cnt - 2);
  dir_vec(:, 1) = 0;
  for i = 1 : keyframe_cnt - 1
    dir_vec(i, i + 1) = 1;
  endfor
  dir_vec_len = sqrt(1 + 1 / (keyframe_cnt - 2));
  dir_vec = dir_vec' / dir_vec_len;

  grad = zeros(keyframe_cnt, 1);
  for i = 1 : keyframe_cnt - 1
    nkl = new_keyframe_list(keyframe_list, dir_vec(:, i), der_len, i);
    [solution, fval] = solve(fourier_order, der_order, nkl, end_point_cond, vel_bound);
    nxt_value = target(fourier_order, solution, fval, keyframe_list, vel_bound);
    grad += dir_vec(:, i) * (cur_value - nxt_value) / der_len;
  endfor
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

function value = target(fourier_order, solution, fval, keyframe_list, vel_bound);
  [min_vel, max_vel] = velocity_range(fourier_order, solution, keyframe_list);
  value = ((max_vel - min_vel) / diff(vel_bound))^6;
  %value = fval;
endfunction
