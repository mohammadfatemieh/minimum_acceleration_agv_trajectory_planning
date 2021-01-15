function [min_value, max_value, avg_speed] = velocity_range(fourier_order, solution, keyframe_list)
  fourier_size = 2 * fourier_order + 1;
  keyframe_cnt = size(keyframe_list, 1);
  time_idx = 3;
  min_value = [];
  max_value = [];
  total_distance = 0;
  for i = 1 : keyframe_cnt - 1
    time_sample = linspace(keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx), 10);
    pos_x_f = fourier(fourier_order, [keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx)]);
    pos_x_f = pos_x_f.scale(solution((i - 1) * fourier_size + 1 : i * fourier_size, 1)');
    vel_x_f = pos_x_f.derivative;
    vel_x = zeros(1, length(time_sample));
    pos_y_f = fourier(fourier_order, [keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx)]);
    pos_y_f = pos_y_f.scale(solution((i - 1 + keyframe_cnt - 1) * fourier_size + 1 : (i + keyframe_cnt - 1) * fourier_size, 1)');
    vel_y_f = pos_y_f.derivative;
    vel_y = zeros(1, length(time_sample));
    for j = 1 : length(time_sample)
      vel_x(j) = sum(vel_x_f.value(time_sample(j)));
      vel_y(j) = sum(vel_y_f.value(time_sample(j)));
      if j > 1
        total_distance += (sqrt(vel_x(j) .^ 2 + vel_y(j) .^ 2) + sqrt(vel_x(j - 1) .^ 2 + vel_y(j - 1) .^ 2)) / 2 * ((keyframe_list(i + 1, time_idx) - keyframe_list(i, time_idx)) / 9);
      endif
    endfor
    min_value = [min_value, min(sqrt(vel_x.^2 + vel_y.^2))];
    max_value = [max_value, max(sqrt(vel_x.^2 + vel_y.^2))];
  endfor
  avg_speed = total_distance / (keyframe_list(keyframe_cnt, time_idx) - keyframe_list(1, time_idx));
  min_value = min(min_value);
  max_value = max(max_value);
endfunction
