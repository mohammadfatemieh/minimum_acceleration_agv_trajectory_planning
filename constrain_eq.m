function [mat, vec] = constrain_eq(poly_order, keyframe_list, vel_cond, acc_cond)

  x_idx = 1;
  y_idx = 2;
  time_idx = size(keyframe_list, 2);
  keyframe_cnt = size(keyframe_list, 1);
  polynomial = ones(1, poly_order);
  vel_poly = polyder(polynomial);
  acc_poly = polyder(vel_poly);
  mat = []; % A
  vec = []; % b

  % For x and y components
  for idx = 1 : 2
    
    %%%%%%%%%%
    %% b
    % path
    pos_vec = keyframe_list(1, idx);
    for i = 2 : keyframe_cnt - 1
      pos_vec = [pos_vec; keyframe_list(i, idx); keyframe_list(i, idx);];
    endfor
    pos_vec = [pos_vec; keyframe_list(keyframe_cnt, idx)];
    % velocity
    vel_vec = [vel_cond(1, idx); zeros((keyframe_cnt - 2), 1); vel_cond(2, idx)];
    % acceleration
    acc_vec = [acc_cond(1, idx); zeros((keyframe_cnt - 2), 1); acc_cond(2, idx)];

    vec = [vec; pos_vec; vel_vec; acc_vec];
    %%%%%%%%%%

    %%%%%%%%%%
    %% A
    L = 2 * (keyframe_list(2, time_idx) - keyframe_list(1, time_idx));
    [pos_A, pos_w, pos_t] = fourier_sin(L, poly_order);
    [vel_A, vel_w, vel_t] = fourier_der(pos_A, pos_w, pos_t);
    [acc_A, acc_w, acc_t] = fourier_der(vel_A, vel_w, vel_t);
    pos_mat = [];
    vel_mat = zeros(1, poly_order * (keyframe_cnt - 1));
    vel_mat(1 : poly_order) = fourier_val(vel_A, vel_w, vel_t, keyframe_list(1, time_idx));
    acc_mat = zeros(1, poly_order * (keyframe_cnt - 1));
    acc_mat(1 : poly_order) = fourier_val(acc_A, acc_w, acc_t, keyframe_list(1, time_idx));
    % position
    for i = 1 : keyframe_cnt - 1
      L = 2 * (keyframe_list(i + 1, time_idx) - keyframe_list(i, time_idx));
      [pos_A, pos_w, pos_t] = fourier_sin(L, poly_order);
      pos_blk = zeros(2, poly_order);
      pos_blk(1, :) = fourier_val(pos_A, pos_w, pos_t, keyframe_list(i, time_idx))';
      pos_blk(2, :) = fourier_val(pos_A, pos_w, pos_t, keyframe_list(i + 1, time_idx))';
      pos_mat = blkdiag(pos_mat, pos_blk);
    endfor
    % velocity
    for i = 1 : keyframe_cnt - 2
      L = 2 * (keyframe_list(i + 1, time_idx) - keyframe_list(i, time_idx));
      [pos_A, pos_w, pos_t] = fourier_sin(L, poly_order);
      [vel_A, vel_w, vel_t] = fourier_der(pos_A, pos_w, pos_t);
      [acc_A, acc_w, acc_t] = fourier_der(vel_A, vel_w, vel_t);
      vel_blk = zeros(1, poly_order * (keyframe_cnt - 1));
      vel_blk((i - 1) * poly_order + 1 : i * poly_order) = fourier_val(vel_A, vel_w, vel_t, keyframe_list(i, time_idx));
      vel_blk(i * poly_order + 1 : (i + 1) * poly_order) = -fourier_val(vel_A, vel_w, vel_t, keyframe_list(i + 1, time_idx));
      vel_blk
      vel_mat = [vel_mat; vel_blk];

      % acceleration
      acc_blk = zeros(1, poly_order * (keyframe_cnt - 1));
      acc_blk((i - 1) * poly_order + 1 : i * poly_order) = fourier_val(acc_A, acc_w, acc_t, keyframe_list(i, time_idx));
      acc_blk(i * poly_order + 1 : (i + 1) * poly_order) = -fourier_val(acc_A, acc_w, acc_t, keyframe_list(i + 1, time_idx));
      acc_blk
      acc_mat = [acc_mat; acc_blk];
    endfor
    vel_mat(size(vel_mat, 1) + 1, (keyframe_cnt - 2) * poly_order + 1 : (keyframe_cnt - 1) * poly_order) = fourier_val(vel_A, vel_w, vel_t, keyframe_list(keyframe_cnt, time_idx));
    acc_mat(size(acc_mat, 1) + 1, (keyframe_cnt - 2) * poly_order + 1 : (keyframe_cnt - 1) * poly_order) = fourier_val(acc_A, acc_w, acc_t, keyframe_list(keyframe_cnt, time_idx));

    sub_mat = [pos_mat; vel_mat; acc_mat];
    mat = blkdiag(mat, sub_mat);
    %%%%%%%%%%

  endfor
endfunction

function substed_poly = poly_sub_t(polynomial, t)
  poly_order = length(polynomial);
  substed_poly = polynomial;
  for i = 1 : poly_order
    substed_poly(i) *= t^(poly_order - i);
  endfor
endfunction
