pkg load optim;

# keyframes = [
#   0, 0, 0;
#   1, 0, 1;
#   1, 2, 4;
#   0, 2, 5;
# ];
keyframes = [
  0, 1, 1;
  0, 4, 4;
];
keyframe_cnt = size(keyframes, 1);

vel_cond = [
  1, 0;
  -1, 0;
];

acc_cond = [
  0, 0;
  0, 0;
];

vel_max = 2;
acc_max = 2;

poly_order = 6;
derivative_order = 2;
cost_blk = zeros(poly_order);

polynomial = ones(1, poly_order);
for i = 1 : derivative_order
  polynomial = polyder(polynomial);
end
cost_blk(1 : poly_order - derivative_order, 1 : poly_order - derivative_order) = polynomial' * polynomial;
for i = 1 : poly_order - derivative_order
  cost_blk(i, 1 : poly_order - i + 1) = polyint(cost_blk(i, 1 : poly_order - i + 1))(1 : poly_order - i + 1);
end

function cost_blk = substitute_t(poly_order, derivative_order, template_cost_blk, cur_t, nxt_t)
  cost_blk = template_cost_blk;
  for i = 1 : poly_order - derivative_order + 1
    for j = 1 : i
      cost_blk(j, i + 1 - j) *= ...
        abs(nxt_t^((poly_order - derivative_order - 1)^2 - i + 2) - cur_t^((poly_order - derivative_order - 1)^2 - i + 2));
    endfor
  endfor
endfunction

time_seq = keyframes(:, 3);
cost_matrix = substitute_t(poly_order, derivative_order, cost_blk, time_seq(1), time_seq(2));
for i = 2 : length(time_seq) - 1
  cost_matrix = blkdiag(cost_matrix, substitute_t(poly_order, derivative_order, cost_blk, time_seq(i), time_seq(i + 1)));
endfor

function [mat, vec] = constrain(poly_order, keyframe_list, vel_cond, acc_cond)
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
    % position
    pos_mat = [];
    for i = 1 : keyframe_cnt - 1
      pos_blk = zeros(2, poly_order);
      pos_blk(1, :) = poly_sub_t(ones(1, poly_order), keyframe_list(i, time_idx));
      pos_blk(2, :) = poly_sub_t(ones(1, poly_order), keyframe_list(i + 1, time_idx));
      pos_mat = blkdiag(pos_mat, pos_blk);
    endfor
    % velocity
    vel_mat = zeros(1, poly_order * (keyframe_cnt - 1));
    vel_mat(1 : poly_order - 1) = poly_sub_t(vel_poly, keyframe_list(1, time_idx));
    for i = 1 : keyframe_cnt - 2
      vel_blk = zeros(1, poly_order * (keyframe_cnt - 1));
      vel_blk((i - 1) * poly_order + 1 : i * poly_order - 1) = poly_sub_t(vel_poly, keyframe_list(i + 1, time_idx));
      vel_blk(i * poly_order + 1 : (i + 1) * poly_order - 1) = -poly_sub_t(vel_poly, keyframe_list(i + 1, time_idx));
      vel_mat = [vel_mat; vel_blk];
    endfor
    vel_mat(size(vel_mat, 1) + 1, (keyframe_cnt - 2) * poly_order + 1 : (keyframe_cnt - 1) * poly_order - 1) = poly_sub_t(vel_poly, keyframe_list(keyframe_cnt, time_idx));
    % acceleration
    acc_mat = zeros(1, poly_order * (keyframe_cnt - 1));
    acc_mat(1 : poly_order - 2) = poly_sub_t(acc_poly, keyframe_list(1, time_idx));
    for i = 1 : keyframe_cnt - 2
      acc_blk = zeros(1, poly_order * (keyframe_cnt - 1));
      acc_blk((i - 1) * poly_order + 1 : i * poly_order - 2) = poly_sub_t(acc_poly, keyframe_list(i + 1, time_idx));
      acc_blk(i * poly_order + 1 : (i + 1) * poly_order - 2) = -poly_sub_t(acc_poly, keyframe_list(i + 1, time_idx));
      acc_mat = [acc_mat; acc_blk];
    endfor
    acc_mat(size(acc_mat, 1) + 1, (keyframe_cnt - 2) * poly_order + 1 : (keyframe_cnt - 1) * poly_order - 2) = poly_sub_t(acc_poly, keyframe_list(keyframe_cnt, time_idx));

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

function sum_array = prefix_sum(array)
  sum_array = array;
  for i = 2 : length(array)
    sum_array(i) += sum_array(i - 1);
  endfor
endfunction

prev_value = 100;
cur_value = 0;

[A, b] = constrain(poly_order, keyframes, vel_cond, acc_cond)
rank(A)
H = blkdiag(cost_matrix, cost_matrix)
f = zeros(2 * (keyframe_cnt - 1) * poly_order, 1);
[result, prev_value] = quadprog(H, f, [], [], A, b)
i = 1
t = linspace(keyframes(i, 3), keyframes(i + 1, 3), 100);
r_x = polyval(result((i - 1) * poly_order + 1 : i * poly_order, 1)', t);
% v_x = polyval(polyder(result((i - 1) * poly_order + 1 : i * poly_order, 1)'), t);
% a_x = polyval(polyder(polyder(result((i - 1) * poly_order + 1 : i * poly_order, 1)')), t);
r_y = polyval(result((i - 1 + keyframe_cnt - 1) * poly_order + 1 : (i + keyframe_cnt - 1) * poly_order, 1)', t);
% v_y = polyval(polyder(result((i - 1 + keyframe_cnt - 1) * poly_order + 1 : (i + keyframe_cnt - 1) * poly_order, 1)'), t);
% a_y = polyval(polyder(polyder(result((i - 1 + keyframe_cnt - 1) * poly_order + 1 : (i + keyframe_cnt - 1) * poly_order, 1)')), t);
plot(r_x, r_y)
% plot(t, (v_x .* a_x + v_y .* a_y) ./ sqrt(v_x.^2 + v_y.^2));
% plot(t, v_x.^2 + v_y.^2)

% while abs(prev_value - cur_value) > 0.01
%   [A, b] = constrain(poly_order, keyframes, vel_cond, acc_cond);
%   H = blkdiag(cost_matrix, cost_matrix);
%   f = zeros(2 * (keyframe_cnt - 1) * poly_order, 1);
%   [result, prev_value] = quadprog(H, f, [], [], A, b);
%   prev_value
%   old_keyframes = keyframes;
%   grad = [0];
%   for i = 1 : keyframe_cnt - 1
%     % step = -ones(keyframe_cnt, 1) / (keyframe_cnt - 2);
%     step = zeros(keyframe_cnt, 1);
%     step(1) = 0;
%     step(i + 1) = 1;
%     keyframes = old_keyframes;
%     keyframes(:, 3) += 0.01 * prefix_sum(step)
%     [A, b] = constrain(poly_order, keyframes, vel_cond, acc_cond);
%     H = blkdiag(cost_matrix, cost_matrix);
%     f = zeros(2 * (keyframe_cnt - 1) * poly_order, 1);
%     [result, cur_value] = quadprog(H, f, [], [], A, b);
%     grad = [grad; (prev_value - cur_value) / 0.01];
%   endfor
%   grad
%   old_keyframes(:, 3) += -0.1 * prefix_sum(grad);
% endwhile
% 
% hold on
% real_roots = [];
% for i = 1 : keyframe_cnt - 1
%   t = linspace(keyframes(i, 3), keyframes(i + 1, 3), 100);
%   r_x = polyval(result((i - 1) * poly_order + 1 : i * poly_order, 1)', t);
%   r_x_poly = result((i - 1) * poly_order + 1 : i * poly_order, 1)';
%   v_x_poly = polyder(r_x_poly);
%   a_x_poly = polyder(v_x_poly);
%   v_x = polyval(polyder(result((i - 1) * poly_order + 1 : i * poly_order, 1)'), t);
%   a_x = polyval(polyder(polyder(result((i - 1) * poly_order + 1 : i * poly_order, 1)')), t);
%   r_y = polyval(result((i - 1 + keyframe_cnt - 1) * poly_order + 1 : (i + keyframe_cnt - 1) * poly_order, 1)', t);
%   r_y_poly = result((i - 1 + keyframe_cnt - 1) * poly_order + 1 : (i + keyframe_cnt - 1) * poly_order, 1)';
%   v_y_poly = polyder(r_y_poly);
%   a_y_poly = polyder(v_y_poly);
%   v_y = polyval(polyder(result((i - 1 + keyframe_cnt - 1) * poly_order + 1 : (i + keyframe_cnt - 1) * poly_order, 1)'), t);
%   a_y = polyval(polyder(polyder(result((i - 1 + keyframe_cnt - 1) * poly_order + 1 : (i + keyframe_cnt - 1) * poly_order, 1)')), t);
%   % plot(r_x, r_y)
%   plot(t, (v_x .* a_x + v_y .* a_y) ./ sqrt(v_x.^2 + v_y.^2));
%   plot(t, v_x.^2 + v_y.^2)
%   poly_roots = roots(conv(v_x_poly, a_x_poly) + conv(v_y_poly, a_y_poly));
%   for j = 1 : length(poly_roots)
%     imag_part = imag(poly_roots(j));
%     real_part = real(poly_roots(j));
%     if ~imag_part && real_part >= keyframes(i, 3) && real_part <= keyframes(i + 1, 3)
%       real_roots = [real_roots, real_part];
%     endif
%   endfor
% endfor
% real_roots
% r_x_roots = polyval(result(1 : 5), real_roots(1))
% r_y_roots = polyval(result(21 : 25), real_roots(1))
% hold off
