geometry = [1, 1]
wheel_radius = 0.05
wheel_distance = 1
vel_bound = [0.5, 1]
% format: [r_x, r_y, theta, v_x, v_y, a_x, a_y]
end_point_cond = [0, 0, 1.57, 0, 1, 0, 0; ...
                  0, 6,    0, 0, 1, 0, 0;];
fourier_order = 5;
fourier_size = 2 * fourier_order + 1;
der_order = 2;
time_seq = [0; 6];
keyframe_list = [end_point_cond(:, 1 : 2), time_seq];
keyframe_list = [keyframe_list(1, :); [-2, 1, 2]; [-2, 3, 4]; keyframe_list(2, :)];
keyframe_cnt = size(keyframe_list, 1);
time_idx = 3;

hold on
while 1
  [coef_list, keyframe_list] = optimize(fourier_order, der_order, keyframe_list, end_point_cond, vel_bound);
  break;
  %min_value = [];
  %max_value = [];
  %for i = 1 : keyframe_cnt - 1
  %  time_sample = linspace(keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx), 50);
  %  pos_x_f = fourier(fourier_order, [keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx)]);
  %  pos_x_f = pos_x_f.scale(coef_list((i - 1) * fourier_size + 1 : i * fourier_size, 1)');
  %  vel_x_f = pos_x_f.derivative;
  %  vel_x = zeros(1, length(time_sample));
  %  pos_y_f = fourier(fourier_order, [keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx)]);
  %  pos_y_f = pos_y_f.scale(coef_list((i - 1 + keyframe_cnt - 1) * fourier_size + 1 : (i + keyframe_cnt - 1) * fourier_size, 1)');
  %  vel_y_f = pos_y_f.derivative;
  %  vel_y = zeros(1, length(time_sample));
  %  for j = 1 : length(time_sample)
  %    vel_x(j) = sum(vel_x_f.value(time_sample(j)));
  %    vel_y(j) = sum(vel_y_f.value(time_sample(j)));
  %  endfor
  %  min_value = [min_value, min(sqrt(vel_x.^2 + vel_y.^2))];
  %  max_value = [max_value, max(sqrt(vel_x.^2 + vel_y.^2))];
  %endfor
  %min_value = min(min_value)
  %max_value = max(max_value)
  %if max_value > 1.2
  %  for j = 1 : keyframe_cnt
  %    keyframe_list(j, time_idx) *= max_value / 1.2;
  %  endfor
  %else
  %  break;
  %endif
  %input('');
endwhile
hold off
