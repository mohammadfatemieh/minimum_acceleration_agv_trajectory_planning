geometry = [1, 1]
wheel_radius = 0.05
wheel_distance = 1
vel_bound = [0.5, 1]
% format: [r_x, r_y, theta, v_x, v_y, a_x, a_y]
end_point_cond = [0, 0, 1.57, 0, 1, 0, 0; ...
                  0, 6,    0, 0, 1, 0, 0;];
poly_order = 6;
der_order = 2;
time_seq = [0; 12];
keyframe_list = [end_point_cond(:, 1 : 2), time_seq];
keyframe_list = [keyframe_list(1, :); [-2, 1, 2]; [-2, 2, 4]; [0, 3, 6]; [2, 4, 8]; [2, 5, 10]; keyframe_list(2, :)];
keyframe_cnt = size(keyframe_list, 1);
time_idx = 3;
low_cnt = 1;
high_cnt = 1;

hold on
while low_cnt ~= 0 || high_cnt ~= 0
  coef_list = solve(poly_order, der_order, keyframe_list, end_point_cond,vel_bound);
  break
  for i = 1 : keyframe_cnt - 1
    t_seq = linspace(keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx), 1000);
    % x
    r_x_poly = coef_list((i - 1) * poly_order + 1 : i * poly_order, 1)';
    v_x_poly = polyder(r_x_poly);
    a_x_poly = polyder(v_x_poly);
    r_x = polyval(r_x_poly, t_seq);
    v_x = polyval(v_x_poly, t_seq);
    a_x = polyval(a_x_poly, t_seq);
    % y
    r_y_poly = coef_list((i - 1 + keyframe_cnt - 1) * poly_order + 1 : (i + keyframe_cnt - 1) * poly_order, 1)';
    v_y_poly = polyder(r_y_poly);
    a_y_poly = polyder(v_y_poly);
    r_y = polyval(r_y_poly, t_seq);
    v_y = polyval(v_y_poly, t_seq);
    a_y = polyval(a_y_poly, t_seq);
    vel = @(t) sqrt(polyval(conv(v_x_poly, v_x_poly) + conv(v_y_poly, v_y_poly), t));
    l = integral(vel, keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx));
    % plot(t_seq, sqrt(v_x.^2 + v_y.^2))
    plot(r_x, r_y)
    real_roots = [];
    poly_roots = roots(conv(v_x_poly, a_x_poly) + conv(v_y_poly, a_y_poly));
    for j = 1 : length(poly_roots)
      imag_part = imag(poly_roots(j));
      real_part = real(poly_roots(j));
      if ~imag_part && real_part >= keyframe_list(i, 3) && real_part <= keyframe_list(i + 1, 3)
        real_roots = [real_roots, real_part];
      endif
    endfor
    real_roots = [keyframe_list(1, time_idx), real_roots, keyframe_list(keyframe_cnt, time_idx)];
  endfor
endwhile
hold off
