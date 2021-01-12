
fourier_size = 2 * fourier_order + 1;
keyframe_cnt = size(keyframe_list, 1);
time_idx = 3;
hold on
real_roots = [];
for i = 1 : keyframe_cnt - 1
  time_sample = linspace(keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx), 50);
  pos_x_f = fourier(fourier_order, [keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx)]);
  pos_x_f = pos_x_f.scale(coef_list((i - 1) * fourier_size + 1 : i * fourier_size, 1)');
  pos_x = zeros(1, length(time_sample));
  for j = 1 : length(time_sample)
    pos_x(j) = sum(pos_x_f.value(time_sample(j)));
  endfor
  vel_x_f = pos_x_f.derivative;
  vel_x = zeros(1, length(time_sample));
  % y
  pos_y_f = fourier(fourier_order, [keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx)]);
  pos_y_f = pos_y_f.scale(coef_list((i - 1 + keyframe_cnt - 1) * fourier_size + 1 : (i + keyframe_cnt - 1) * fourier_size, 1)');
  pos_y = zeros(1, length(time_sample));
  for j = 1 : length(time_sample)
    pos_y(j) = sum(pos_y_f.value(time_sample(j)));
  endfor
  vel_y_f = pos_y_f.derivative;
  vel_y = zeros(1, length(time_sample));
  for j = 1 : length(time_sample)
    vel_x(j) = sum(vel_x_f.value(time_sample(j)));
    vel_y(j) = sum(vel_y_f.value(time_sample(j)));
  endfor
  plot(time_sample, sqrt(vel_x.^2 + vel_y.^2));
  % plot(vel_x_f);
  % input('');
  % plot(vel_y_f);
  % plot(time_sample, v_x);
  % plot(time_sample, v_y);
  % plot(r_x, r_y);
endfor
hold off
