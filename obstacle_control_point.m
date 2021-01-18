function control_point_list  = obstacle_control_point(obstacle, safe_distance, start_point)
  % control_point_list  Generate control point list from obstacle, safe_distance and start/end point.
  obstacle_vec = [(obstacle(2, 1) - obstacle(1, 1)), (obstacle(2, 2) - obstacle(1, 2))];
  start_obstacle_vec = [(start_point(1) - obstacle(1, 1)), (start_point(2) - obstacle(1, 2))];
  start_cross_obstacle = cross2d(start_obstacle_vec, obstacle_vec);
  obstacle_complex = obstacle_vec(1) + obstacle_vec(2) * 1i;
  obstacle_complex *= 1i;
  new_obstacle_vec = [real(obstacle_complex), imag(obstacle_complex)];
  obstacle_cross_obstacle = cross2d(new_obstacle_vec, obstacle_vec);
  if obstacle_cross_obstacle * start_cross_obstacle > 0
    new_obstacle_vec = -new_obstacle_vec;
  endif
  dangerous_region = [obstacle; obstacle(1, :) + new_obstacle_vec; obstacle(2, :) + new_obstacle_vec];
  dangerous_region_center = sum(dangerous_region) / 4;
  for i = 4 : -1 : 1
    center_vertex_vec = (dangerous_region(i, :) - dangerous_region_center);
    center_vertex_vec = center_vertex_vec / norm(center_vertex_vec);
    dangerous_region(i, :) = dangerous_region(mod(i + 1, 2) + 1, :) + safe_distance * center_vertex_vec;
    figure(1);
    hold on;
    plot(dangerous_region(i, 1), dangerous_region(i, 2), 'ro');
    hold off;
  endfor
  control_point_list = dangerous_region;
endfunction

function product = cross2d(vec1, vec2)
  product = vec1(1) * vec2(2) - vec1(2) * vec2(1);
endfunction
