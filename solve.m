function [solution, value] = solve(fourier_order, der_order, keyframe_list, end_point_cond, vel_bound)
  fourier_size = 2 * fourier_order + 1;
  keyframe_cnt = size(keyframe_list, 1);
  time_idx = 3;
  vel_cond = end_point_cond(1 : 2, 4 : 5);
  acc_cond = end_point_cond(1 : 2, 6 : 7);

  H = cost(fourier_order, der_order, keyframe_list);
  f = zeros(2 * (keyframe_cnt - 1) * fourier_size, 1);
  [A, b] = constrain_eq(fourier_order, keyframe_list, vel_cond, acc_cond);
  %[Aieq, bieq] = constrain_ieq(fourier_order, keyframe_list);
  %[solution, value] = quadprog(H, f, [], [], A, b, -3 * ones(2 * (keyframe_cnt - 1) * fourier_size, 1), 3 * ones(2 * (keyframe_cnt - 1) * fourier_size, 1));
  [solution, value] = quadprog(H, f, [], [], A, b);
  %[solution, value] = quadprog(H, f, Aieq, bieq, A, b);
endfunction
