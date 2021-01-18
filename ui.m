% - **曲线路径规划**：给定一个双轮差速驱动机器人，和路径的起止点，实现机器人的曲线路径规划，包括轨迹规划、速度规划等；
%   - 输入：
%     - 机器人参数，如大小、轮半径、左右轮的轮距，最大/最小速度；
%     - 路径起止点 (x, y) 和/或障碍物信息（障碍物的横向大小）：
%     - 机器人当前位姿 (x, y, theta)，和左右轮速度；
%   - 输出：
%     - 机器人当前时刻的速度（矢量）；
%   - 要求：
%     - 根据输入信息，完成机器人的曲线路径规划，要求曲线（分段）光滑，一阶可微（即速度连续）；
%     - 根据输入信息，完成机器人的加速度、速度规划，给出当前时刻机器人应有的速度；
%     - **改变目标点，实时进行机器人的路径重规划，并给出当前时刻的速度；**
%     - 实现语言： c++, python, matlab 均可；
%   - 参考文献：
%     - Minimum Snap Trajectory Generation and Control for Quadrotors

global vel_bound end_point_vel end_point_theta;

%%% UI 相关内容
%% 读取用户输入
geometry = [1, 1];
geometry = input('geometry (e.g. [x, y]): ');
vel_bound = [0.8, 1.2];
vel_bound = input('maximum velocity: ');
vel_bound = [0, vel_bound];
start_point_theta = input('theta at the start point (in degree): ');
start_point_theta = start_point_theta / 180 * pi;
end_point_theta = input('theta at the end point (in degree): ');
end_point_theta = end_point_theta / 180 * pi;
start_point_vel = input('velocity at the start point (magnitude): ');
start_point_vel = max(0.1, start_point_vel);
end_point_vel = input('velocity at the end point (magnitude): ');
end_point_vel = max(0.1, end_point_vel);
start_point_vel = [cos(start_point_theta) * start_point_vel, sin(start_point_theta) * start_point_vel];
end_point_vel = [cos(end_point_theta) * end_point_vel, sin(end_point_theta) * end_point_vel];

figure(1);

% 设置障碍的回调函数
function set_obstacle_cb(h, ~)
  global is_obstacle obstacle obstacle_done;
  if obstacle_done
    return;
  endif
  pos = get(gca, 'CurrentPoint');
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

% 设定起点和终点的回调函数
function set_start_end_cb(h, ~)
  global start_or_end start_point end_point start_end_done;
  if start_end_done
    return;
  endif
  pos = get(gca, 'CurrentPoint');
  if start_or_end == logical(0) % start
    start_point = pos(1, 1 : 2);
  else
    end_point = pos(1, 1 : 2);
    start_end_done = logical(1);
  endif
  start_or_end = ~start_or_end;
  plot(pos(1, 1), pos(1, 2), 'bo');
endfunction

% 初始化用户界面
grid on; hold on;
xlabel('X-Y Plane');
axis([-10 10 -10 10]);
global is_obstacle obstacle obstacle_done;
obstacle = zeros(2);
is_obstacle = logical(0);
obstacle_done = logical(0);

% 设置障碍
set(gcf, 'WindowButtonDownFcn', @set_obstacle_cb);

% 等待用户完成
disp('Start obstacle setting.');
while obstacle_done ~= logical(1)
  drawnow;
endwhile
disp('Done obstacle setting.');
obstacle_vec = [(obstacle(2, 1) - obstacle(1, 1)), (obstacle(2, 2) - obstacle(1, 2))];

% 设定起止点
global start_or_end start_point end_point start_end_done;
start_or_end = logical(0);
start_point = zeros(1, 2);
end_point = zeros(1, 2);
start_end_done = logical(0);

set(gcf, 'WindowButtonDownFcn', @set_start_end_cb);

% 等待用户完成
disp('Start start and end point setting.');
while start_end_done ~= logical(1)
  drawnow;
endwhile
disp('Done start and end point setting.');
disp('Start planning...');

drawnow;

global keyframe_list;

% 初始化参数
fourier_order = 3;
fourier_size = 2 * fourier_order + 1;
% 最小化目标相对位移的求道阶数
% 即:
%   1 -> minimum speed
%   2 -> minimum acceleration
%   3 -> minimum jerk
%   4 -> minimum snap
der_order = 2;
time_idx = 3;
start_point = [start_point(1, :), start_point_theta, start_point_vel, 0, 0];
end_point = [end_point(1, :), end_point_theta, end_point_vel, 0, 0];
assumed_velocity = sum(vel_bound) / 2;

% 设定安全距离
safe_distance = norm(geometry) * 1.5;

global dangerous_region;

% 主循环
for cnt = 0 : 2
  control_point_list = obstacle_control_point(obstacle, safe_distance, start_point);
  dangerous_region = obstacle_control_point(obstacle, norm(geometry), start_point);
  [keyframe_list, end_point_cond] = init_plan( ...
    control_point_list, ...
    cnt, ...
    start_point, ...
    end_point, ...
    vel_bound, ...
    assumed_velocity);
  [solution, _] = solve(fourier_order, der_order, keyframe_list, end_point_cond, vel_bound);
  [min_vel, max_vel, avg_speed] = velocity_range(fourier_order, solution, keyframe_list)
  if verify(fourier_order, solution, keyframe_list, dangerous_region)
    if max_vel > vel_bound(2)
      keyframe_list(:, 3) *= max_vel / vel_bound(2);
      for k = 1 : 5
        [solution, _] = solve(fourier_order, der_order, keyframe_list, end_point_cond, vel_bound);
        [min_vel, max_vel, avg_speed] = velocity_range(fourier_order, solution, keyframe_list)
        if max_vel > vel_bound(2)
          keyframe_list(:, 3) *= max_vel / vel_bound(2) * 1.05;
        else
          break;
        endif
      endfor
      if max_vel > vel_bound(2)
        [solution, keyframe_list] = optimize(fourier_order, der_order, keyframe_list, end_point_cond, vel_bound);
        [min_vel, max_vel, avg_speed] = velocity_range(fourier_order, solution, keyframe_list);
        for k = 1 : 5
          if max_vel > vel_bound(2)
            keyframe_list(:, 3) *= max_vel / vel_bound(2) * 1.05;
          else
            break;
          endif
          [solution, _] = solve(fourier_order, der_order, keyframe_list, end_point_cond, vel_bound);
          [min_vel, max_vel, avg_speed] = velocity_range(fourier_order, solution, keyframe_list);
        endfor
      else
        if verify(fourier_order, solution, keyframe_list, dangerous_region)
          break;
        else
          continue;
        endif
      endif
    else
      break;
    endif
  endif
endfor
plot_traj(fourier_order, solution, keyframe_list);

global current_time current_pos current_vel current_stage;

% 改变终点实时规划新路径
% boilerplot of main loop {{{
function change_end_point_cb(h, ~)
  global start_point end_point obstacle;
  global current_time current_pos current_vel current_stage;
  global keyframe_list vel_bound end_point_vel end_point_theta;
  global vel_x vel_y time_sample;
  global dangerous_region;
  figure(2);
  hold on;
  plot(time_sample, sqrt(vel_x .^ 2 + vel_y .^ 2));
  drawnow;
  hold off;
  click_type = get(h, 'selectiontype');
  figure(1);
  pos = get(gca, 'CurrentPoint');
  if strcmp(click_type, 'normal')
    figure(1);
    hold on;
    plot(pos(1, 1), pos(1, 2), 'bo');
    drawnow;
    hold off;
    fourier_order = 5;
    fourier_size = 2 * fourier_order + 1;
    der_order = 2;
    time_idx = 3;
    assumed_velocity = sum(vel_bound) / 2
    safe_distance = 1;
    start_point = [current_pos, 0, current_vel, 0, 0]
    end_point = [pos(1, 1 : 2), end_point_theta, end_point_vel, 0, 0]
    end_point_cond = [start_point; end_point];
    old_keyframe_list = keyframe_list;
    for cnt = 0 : size(old_keyframe_list, 1) - 1 - current_stage
      keyframe_list = zeros(cnt + 2, 3);
      keyframe_list(:, 1 : 2) = [start_point(1 : 2); old_keyframe_list(current_stage + 1 : current_stage + cnt, 1 : 2); end_point(1 : 2);];
      keyframe_list(1, 3) = 0;
      for i = 2 : size(keyframe_list, 1)
        keyframe_list(i, 3) = sqrt(sum((keyframe_list(i, 1 : 2) - keyframe_list(i - 1, 1 : 2)).^2)) / assumed_velocity + keyframe_list(i - 1, 3);
      endfor
      [solution, _] = solve(fourier_order, der_order, keyframe_list, end_point_cond, vel_bound);
      [min_vel, max_vel, avg_speed] = velocity_range(fourier_order, solution, keyframe_list);
      if verify(fourier_order, solution, keyframe_list, dangerous_region)
        if max_vel > vel_bound(2)
          keyframe_list(:, 3) *= max_vel / vel_bound(2);
          for k = 1 : 3
            [solution, _] = solve(fourier_order, der_order, keyframe_list, end_point_cond, vel_bound);
            [min_vel, max_vel, avg_speed] = velocity_range(fourier_order, solution, keyframe_list)
            if max_vel > vel_bound(2)
              keyframe_list(:, 3) *= max_vel / vel_bound(2);
            else
              break;
            endif
          endfor
          if max_vel > vel_bound(2)
            [solution, keyframe_list] = optimize(fourier_order, der_order, keyframe_list, end_point_cond, vel_bound);
            [min_vel, max_vel, avg_speed] = velocity_range(fourier_order, solution, keyframe_list);
            for k = 1 : 3
              if max_vel > vel_bound(2)
                keyframe_list(:, 3) *= max_vel / vel_bound(2);
              else
                break;
              endif
              [solution, _] = solve(fourier_order, der_order, keyframe_list, end_point_cond, vel_bound);
              [min_vel, max_vel, avg_speed] = velocity_range(fourier_order, solution, keyframe_list);
            endfor
          else
            if verify(fourier_order, solution, keyframe_list, dangerous_region)
              break;
            else
              continue;
            endif
          endif
        else
          break;
        endif
      endif
    endfor
    keyframe_list(:, 3) += current_time;
    figure(1);
    hold on;
    plot_traj(fourier_order, solution, keyframe_list);
    hold off;
    keyframe_cnt = size(keyframe_list, 1);
    for i = 1 : keyframe_cnt - 1
      current_stage = i;
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
        pos_x = sum(pos_x_f.value(time_sample(j)));
        pos_y = sum(pos_y_f.value(time_sample(j)));
        current_time = time_sample(j);
        current_pos = [pos_x, pos_y];
        current_vel = [vel_x(j), vel_y(j)];
        disp([current_vel, norm(current_vel)]);
        figure(1);
        hold on;
        plot(sum(pos_x_f.value(time_sample(j))), sum(pos_y_f.value(time_sample(j))), 'bo');
        drawnow;
        hold off;
        input('');
      endfor
      figure(2);
      hold on;
      plot(time_sample, sqrt(vel_x .^ 2 + vel_y .^ 2));
      hold off;
    endfor
    hold off;
    while logical(1)
      input('');
    end
    exit;
  endif
endfunction

set(gcf, 'WindowButtonDownFcn', @change_end_point_cb);

keyframe_cnt = size(keyframe_list, 1);
hold on;
for i = 1 : keyframe_cnt - 1
  global time_sample vel_x vel_y;
  current_stage = i;
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
    pos_x = sum(pos_x_f.value(time_sample(j)));
    pos_y = sum(pos_y_f.value(time_sample(j)));
    current_time = time_sample(j);
    current_pos = [pos_x, pos_y];
    current_vel = [vel_x(j), vel_y(j)];
    disp([current_vel, norm(current_vel)]);
    figure(1);
    hold on;
    plot(sum(pos_x_f.value(time_sample(j))), sum(pos_y_f.value(time_sample(j))), 'bo');
    drawnow;
    hold off;
    input('');
  endfor
  figure(2);
  hold on;
  plot(time_sample, sqrt(vel_x .^ 2 + vel_y .^ 2));
  hold off;
endfor
hold off;
% }}}
