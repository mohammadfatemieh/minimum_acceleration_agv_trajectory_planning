function [mat, vec] = constrain_eq(fourier_order, keyframe_list, vel_cond, acc_cond)
  % constrain_eq  Return the equality constrain in the form of [A, b]

  x_idx = 1;
  y_idx = 2;
  time_idx = size(keyframe_list, 2);
  keyframe_cnt = size(keyframe_list, 1);
  mat = []; % A
  vec = []; % b

  % For both x and y components
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
    fourier_size = 2 * fourier_order + 1;
    pos_f = fourier(fourier_order, [keyframe_list(1, time_idx), keyframe_list(2, time_idx)]);
    vel_f = pos_f.derivative;
    acc_f = vel_f.derivative;
    pos_mat = [];
    % position
    % At every keyframe, we need to let the position of the AGV to be exactly the same as the designated one.
    % Therefore, we need A and b to be like
    % A = 
    %   [ pos_f_1(t0) ], 0, ...
    %   [ pos_f_1(t1) ], 0, ...
    %   0, ...        0, [ pos_f_2(t1) ], 0, ...
    %   0, ...        0, [ pos_f_2(t2) ], 0, ...
    %   ...
    % b = [x(t0), x(t1), x(t1), x(t2), ...]
    for i = 1 : keyframe_cnt - 1
      pos_f = fourier(fourier_order, [keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx)]);
      pos_blk = zeros(2, fourier_size);
      for j = 1 : 2
        pos_blk(j, :) = pos_f.value(keyframe_list(i + j - 1, time_idx));
      endfor
      pos_mat = blkdiag(pos_mat, pos_blk);
    endfor
    % velocity
    % We want to have smooth transition at the keyframe point, therefore for each row
    % ..., 0, [ vel_f_1(t1) ], [ -vel_f_2(t1) ], 0, ...
    % and corresponding row in b should be 0
    vel_mat = zeros(1, fourier_size * (keyframe_cnt - 1));
    vel_mat(1 : fourier_size) = vel_f.value(keyframe_list(1, time_idx));
    acc_mat = zeros(1, fourier_size * (keyframe_cnt - 1));
    acc_mat(1 : fourier_size) = acc_f.value(keyframe_list(1, time_idx));
    for i = 2 : keyframe_cnt - 1
      pos_f = {fourier(fourier_order, [keyframe_list(i - 1, time_idx), keyframe_list(i - 0, time_idx)]), ...
               fourier(fourier_order, [keyframe_list(i + 0, time_idx), keyframe_list(i + 1, time_idx)])};
      vel_f = {pos_f{1}.derivative, pos_f{2}.derivative};
      acc_f = {vel_f{1}.derivative, vel_f{2}.derivative};
      vel_blk = zeros(1, fourier_size * (keyframe_cnt - 1));
      vel_blk((i - 2) * fourier_size + 1 : (i - 1) * fourier_size) =  vel_f{1}.value(keyframe_list(i, time_idx));
      vel_blk((i - 1) * fourier_size + 1 : (i + 0) * fourier_size) = -vel_f{2}.value(keyframe_list(i, time_idx));
      vel_mat = [vel_mat; vel_blk];

      % acceleration
      acc_blk = zeros(1, fourier_size * (keyframe_cnt - 1));
      acc_blk((i - 2) * fourier_size + 1 : (i - 1) * fourier_size) =  acc_f{1}.value(keyframe_list(i, time_idx));
      acc_blk((i - 1) * fourier_size + 1 : (i + 0) * fourier_size) = -acc_f{2}.value(keyframe_list(i, time_idx));
      acc_mat = [acc_mat; acc_blk];
    endfor
    pos_f = fourier(fourier_order, [keyframe_list(keyframe_cnt - 1, time_idx), keyframe_list(keyframe_cnt, time_idx)]);
    vel_f = pos_f.derivative;
    acc_f = vel_f.derivative;
    vel_mat(size(vel_mat, 1) + 1, (keyframe_cnt - 2) * fourier_size + 1 : (keyframe_cnt - 1) * fourier_size) = vel_f.value(keyframe_list(keyframe_cnt, time_idx));
    acc_mat(size(acc_mat, 1) + 1, (keyframe_cnt - 2) * fourier_size + 1 : (keyframe_cnt - 1) * fourier_size) = acc_f.value(keyframe_list(keyframe_cnt, time_idx));

    sub_mat = [pos_mat; vel_mat; acc_mat];
    mat = blkdiag(mat, sub_mat);
    %%%%%%%%%%

  endfor
endfunction
