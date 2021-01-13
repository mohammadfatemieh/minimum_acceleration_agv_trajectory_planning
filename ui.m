pkg load optim;
figure(1);

function set_obstacle_cb(h, ~)
  global is_obstacle obstacle obstacle_done;
  if obstacle_done
    return;
  endif
  pos = get(gca, 'CurrentPoint')
  if is_obstacle
    obstacle(2, :) = pos(1, 1 : 2);
    plot(obstacle(:, 1), obstacle(:, 2));
    obstacle_done = logical(1);
  else
    obstacle(1, :) = pos(1, 1 : 2);
  endif
  is_obstacle = ~is_obstacle;
  plot(pos(1, 1), pos(1, 2), 'ro');
endfunction

function set_start_end_cb(h, ~)
  global start_or_end start_point end_point start_end_done;
  if start_end_done
    return;
  endif
  pos = get(gca, 'CurrentPoint')
  if start_or_end == logical(0) % start
    start_point = pos(1, 1 : 2);
  else
    end_point = pos(1, 1 : 2);
    start_end_done = logical(1);
  endif
  start_or_end = ~start_or_end;
  plot(pos(1, 1), pos(1, 2), 'bo');
endfunction

function product = cross2d(vec1, vec2)
  product = vec1(1) * vec2(2) - vec1(2) * vec2(1);
endfunction

function result = is_intersected(seg_a, seg_b)
  if min(seg_a(:, 1)) <= max(seg_b(:, 1)) && ...
     min(seg_b(:, 1)) <= max(seg_a(:, 1)) && ...
     min(seg_a(:, 2)) <= max(seg_b(:, 2)) && ...
     min(seg_b(:, 2)) <= max(seg_a(:, 2)) && ...
     (cross2d(seg_a(1, :) - seg_b(1, :), seg_b(2, :) - seg_b(1, :)) * ...
      cross2d(seg_a(2, :) - seg_b(1, :), seg_b(2, :) - seg_b(1, :)) <= 0)
    result = logical(1);
  else
    result = logical(0);
  endif
endfunction

grid on; hold on;
xlabel('X-Y Plane');
axis([-5 20 -12.5 12.5]);
global is_obstacle obstacle obstacle_done;
obstacle = zeros(2);
is_obstacle = logical(0);
obstacle_done = logical(0);
set(gcf, 'WindowButtonDownFcn', @set_obstacle_cb);

while obstacle_done ~= logical(1)
  drawnow;
endwhile
disp('Done obstacle setting.');
obstacle_vec = [(obstacle(2, 1) - obstacle(1, 1)), (obstacle(2, 2) - obstacle(1, 2))];

global start_or_end start_point end_point start_end_done;
start_or_end = logical(0);
start_point = zeros(1, 2);
end_point = zeros(1, 2);
start_end_done = logical(0);
set(gcf, 'WindowButtonDownFcn', @set_start_end_cb);

while start_end_done ~= logical(1)
  drawnow;
endwhile
disp('Done start and end point setting.');

% start_obstacle_vec = [(start_point(1) - obstacle(1, 1)), (start_point(2) - obstacle(1, 2))];
% start_cross_obstacle = cross2d(start_obstacle_vec, obstacle_vec);
% obstacle_complex = obstacle_vec(1) + obstacle_vec(2) * 1i;
% obstacle_complex *= 1i;
% new_obstacle_vec = [real(obstacle_complex), imag(obstacle_complex)]
% obstacle_cross_obstacle = cross2d(new_obstacle_vec, obstacle_vec);
% if obstacle_cross_obstacle * start_cross_obstacle > 0
%   new_obstacle_vec = -new_obstacle_vec;
% endif
% dangerous_region = [obstacle; obstacle(1, :) + new_obstacle_vec; obstacle(2, :) + new_obstacle_vec];
% dangerous_region_center = sum(dangerous_region) / 4;
% for i = 4 : -1 : 1
%   center_vertex_vec = (dangerous_region(i, :) - dangerous_region_center);
%   center_vertex_vec = center_vertex_vec / norm(center_vertex_vec);
%   dangerous_region(i, :) = dangerous_region(mod(i + 1, 2) + 1, :) + 2 * center_vertex_vec;
%   %dangerous_region(i, :) += center_vertex_vec;
%   plot(dangerous_region(i, 1), dangerous_region(i, 2), 'ro');
% endfor

drawnow;

fourier_order = 5;
fourier_size = 2 * fourier_order + 1;
der_order = 2;
time_idx = 3;
start_point = [start_point(1, :), 0, 0, 1, 0, 0];
end_point = [end_point(1, :),   0, 0, 1, 0, 0];
vel_bound = [0.8, 1];
assumed_velocity = sum(vel_bound) / 2;

safe_distance = 1;

for k = 1 : 100
  [keyframe_list, end_point_cond] = init_plan(obstacle, safe_distance, start_point, end_point, vel_bound, assumed_velocity);
  [solution, keyframe_list] = optimize(fourier_order, der_order, keyframe_list, end_point_cond, vel_bound);
  plot_traj(fourier_order, solution, keyframe_list);
  [min_vel, max_vel] = velocity_range(fourier_order, solution, keyframe_list)
  if (max_vel - min_vel) / (max_vel + min_vel) < diff(vel_bound) / sum(vel_bound)
    break;
  else
    safe_distance += 0.2;
  endif
endfor

%end_point_cond = [start_point(1, :), 0, 0, 1, 0, 0;
%                  end_point(1, :),   0, 0, 1, 0, 0;];
%vel_bound = [0.8, 1.2];
%assumed_velocity = sum(vel_bound) / length(vel_bound);
%keyframe_list = [start_point, 0; end_point, norm(diff(end_point)) / assumed_velocity]
%keyframe_cnt = size(keyframe_list, 1);
%[solution, keyframe_list] = optimize(fourier_order, der_order, keyframe_list, end_point_cond, [0.5, 1.5]);
%plot_traj(fourier_order, solution, keyframe_list);
%drawnow;
%
%if verify(fourier_order, solution, keyframe_list, dangerous_region) == 0
%  keyframe_list = zeros(4, 3, 2);
%
%  for k = 1 : 2
%    obstacle_control_point = [dangerous_region(k, :);
%                              dangerous_region(k + 2, :)];
%    keyframe_list(:, 1 : 2, k) = [start_point; obstacle_control_point; end_point];
%
%    % Rough plan of time
%    keyframe_list(1, 3, k) = 0;
%    for i = 2 : size(keyframe_list, 1)
%      keyframe_list(i, 3, k) = sqrt(sum((keyframe_list(i, 1 : 2, k) - keyframe_list(i - 1, 1 : 2, k)).^2)) / assumed_velocity + keyframe_list(i - 1, 3, k);
%    endfor
%  endfor
%
%  if keyframe_list(4, 3, 1) < keyframe_list(4, 3, 2)
%    keyframe_list = keyframe_list(:, :, 1);
%  else
%    keyframe_list = keyframe_list(:, :, 2);
%  endif
%
%  [solution, _] = solve(fourier_order, der_order, keyframe_list, end_point_cond, vel_bound);
%  [min_vel, max_vel] = velocity_range(fourier_order, solution, keyframe_list);
%  if max_vel - min_vel > diff(vel_bound)
%  endif
%
%  [solution, keyframe_list] = optimize(fourier_order, der_order, keyframe_list, end_point_cond, [0.5, 1.5]);
%  plot_traj(fourier_order, solution, keyframe_list);
%  drawnow;
%endif

keyframe_cnt = size(keyframe_list, 1);
hold on;
for i = 1 : keyframe_cnt - 1
  time_sample = linspace(keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx), round((keyframe_list(i + 1, time_idx) - keyframe_list(i, time_idx)) / 0.1));
  pos_x_f = fourier(fourier_order, [keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx)]);
  pos_x_f = pos_x_f.scale(solution((i - 1) * fourier_size + 1 : i * fourier_size, 1)');
  pos_y_f = fourier(fourier_order, [keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx)]);
  pos_y_f = pos_y_f.scale(solution((i - 1 + keyframe_cnt - 1) * fourier_size + 1 : (i + keyframe_cnt - 1) * fourier_size, 1)');
  vel_x_f = pos_x_f.derivative;
  vel_x = zeros(1, length(time_sample));
  vel_y_f = pos_y_f.derivative;
  vel_y = zeros(1, length(time_sample));
  for j = 1 : length(time_sample)
    vel_x(j) = sum(vel_x_f.value(time_sample(j)));
    vel_y(j) = sum(vel_y_f.value(time_sample(j)));
    pos_x = sum(pos_x_f.value(time_sample(j)))
    pos_y = sum(pos_y_f.value(time_sample(j)))
    figure(1);
    hold on;
    plot(sum(pos_x_f.value(time_sample(j))), sum(pos_y_f.value(time_sample(j))), 'bo');
    drawnow;
    hold off;
    input('');
  endfor
  figure(2);
  subplot(2, 1, 1);
  hold on;
  plot(time_sample, vel_x);
  hold off;
  subplot(2, 1, 2);
  hold on;
  plot(time_sample, vel_y);
  hold off;
  figure(3);
  hold on;
  plot(time_sample, sqrt(vel_x .^ 2 + vel_y .^ 2));
  hold off;
endfor
hold off;
