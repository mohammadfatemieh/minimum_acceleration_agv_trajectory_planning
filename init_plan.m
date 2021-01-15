function [keyframe_list_with_time, end_point_cond, dangerous_region] = init_plan(control_point_list, control_point_cnt, start_point, end_point, vel_bound, assumed_velocity)
  %obstacle_vec = [(obstacle(2, 1) - obstacle(1, 1)), (obstacle(2, 2) - obstacle(1, 2))];
  %start_obstacle_vec = [(start_point(1) - obstacle(1, 1)), (start_point(2) - obstacle(1, 2))];
  %start_cross_obstacle = cross2d(start_obstacle_vec, obstacle_vec);
  %obstacle_complex = obstacle_vec(1) + obstacle_vec(2) * 1i;
  %obstacle_complex *= 1i;
  %new_obstacle_vec = [real(obstacle_complex), imag(obstacle_complex)]
  %obstacle_cross_obstacle = cross2d(new_obstacle_vec, obstacle_vec);
  %if obstacle_cross_obstacle * start_cross_obstacle > 0
  %  new_obstacle_vec = -new_obstacle_vec;
  %endif
  %dangerous_region = [obstacle; obstacle(1, :) + new_obstacle_vec; obstacle(2, :) + new_obstacle_vec];
  %dangerous_region_center = sum(dangerous_region) / 4;
  %for i = 4 : -1 : 1
  %  center_vertex_vec = (dangerous_region(i, :) - dangerous_region_center);
  %  center_vertex_vec = center_vertex_vec / norm(center_vertex_vec);
  %  dangerous_region(i, :) = dangerous_region(mod(i + 1, 2) + 1, :) + safe_distance * center_vertex_vec;
  %  figure(1);
  %  hold on;
  %  plot(dangerous_region(i, 1), dangerous_region(i, 2), 'ro');
  %  hold off;
  %endfor
  %drawnow;
  dangerous_region = control_point_list;
  end_point_cond = [start_point(1, :);
                    end_point(1, :)];
  if control_point_cnt == 0
    keyframe_list = zeros(2, 3);
    keyframe_list(:, 1 : 2) = [start_point(1, 1 : 2); end_point(1, 1 : 2)];
    keyframe_list(1, 3) = 0;
    for i = 2 : size(keyframe_list, 1)
      keyframe_list(i, 3) = sqrt(sum((keyframe_list(i, 1 : 2) - keyframe_list(i - 1)).^2)) / assumed_velocity + keyframe_list(i - 1, 3);
    endfor
    keyframe_list_with_time = keyframe_list;
  else
    if control_point_cnt == 1
      keyframe_list = zeros(3, 3, 2);
      for k = 1 : 2
        obstacle_control_point = dangerous_region(k, :);
        keyframe_list(:, 1 : 2, k) = [start_point(1, 1 : 2); obstacle_control_point; end_point(1, 1 : 2)];
        keyframe_list(1, 3, k) = 0;
        for i = 2 : size(keyframe_list, 1)
          keyframe_list(i, 3, k) = sqrt(sum((keyframe_list(i, 1 : 2, k) - keyframe_list(i - 1, 1 : 2, k)).^2)) / assumed_velocity + keyframe_list(i - 1, 3, k);
        endfor
      endfor
      if keyframe_list(3, 3, 1) < keyframe_list(3, 3, 2)
        keyframe_list = keyframe_list(:, :, 1);
      else
        keyframe_list = keyframe_list(:, :, 2);
      endif
      keyframe_list_with_time = keyframe_list;
    else
      keyframe_list = zeros(4, 3, 2);
      for k = 1 : 2
        obstacle_control_point = [dangerous_region(k, :);
                                  dangerous_region(k + 2, :)];
        keyframe_list(:, 1 : 2, k) = [start_point(1, 1 : 2); obstacle_control_point; end_point(1, 1 : 2)];
        keyframe_list(1, 3, k) = 0;
        for i = 2 : size(keyframe_list, 1)
          keyframe_list(i, 3, k) = sqrt(sum((keyframe_list(i, 1 : 2, k) - keyframe_list(i - 1, 1 : 2, k)).^2)) / assumed_velocity + keyframe_list(i - 1, 3, k);
        endfor
      endfor
      if keyframe_list(4, 3, 1) < keyframe_list(4, 3, 2)
        keyframe_list = keyframe_list(:, :, 1);
      else
        keyframe_list = keyframe_list(:, :, 2);
      endif
      keyframe_list_with_time = keyframe_list;
    endif
  endif
endfunction
