function plot_traj(fourier_order, coef_list, keyframe_list)
  % plot_traj  Plot the trajectory
  fourier_size = 2 * fourier_order + 1;
  keyframe_cnt = size(keyframe_list, 1);
  time_idx = 3;
  real_roots = [];
  for i = 1 : keyframe_cnt - 1
    t_seq = linspace(keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx), 50);
    % x
    r_x_f = fourier(fourier_order, [keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx)]);
    r_x_f = r_x_f.scale(coef_list((i - 1) * fourier_size + 1 : i * fourier_size, 1)');
    r_x = zeros(1, length(t_seq));
    for j = 1 : length(t_seq)
      r_x(j) = sum(r_x_f.value(t_seq(j)));
    endfor
    v_x_f = r_x_f.derivative;
    v_x = zeros(1, length(t_seq));
    % for j = 1 : length(t_seq)
    % endfor
    % y
    r_y_f = fourier(fourier_order, [keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx)]);
    r_y_f = r_y_f.scale(coef_list((i - 1 + keyframe_cnt - 1) * fourier_size + 1 : (i + keyframe_cnt - 1) * fourier_size, 1)');
    r_y = zeros(1, length(t_seq));
    for j = 1 : length(t_seq)
      r_y(j) = sum(r_y_f.value(t_seq(j)));
    endfor
    v_y_f = r_y_f.derivative;
    v_y = zeros(1, length(t_seq));
    for j = 1 : length(t_seq)
      v_x(j) = sum(v_x_f.value(t_seq(j)));
      v_y(j) = sum(v_y_f.value(t_seq(j)));
    endfor
    % plot(v_x_f);
    % input('');
    % plot(v_y_f);
    % plot(t_seq, v_x);
    % plot(t_seq, v_y);
    hold on;
    plot(r_x, r_y);
    %plot(t_seq, sqrt(v_x.^2 + v_y.^2));
    hold off;
  endfor
endfunction
