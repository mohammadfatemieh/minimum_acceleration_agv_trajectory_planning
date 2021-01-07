function plot_traj(poly_order, coef_list, keyframe_list)
  keyframe_cnt = size(keyframe_list, 1);
  time_idx = 3;
  hold on
  real_roots = [];
  for i = 1 : keyframe_cnt - 1
    t_seq = linspace(keyframe_list(i, time_idx), keyframe_list(i + 1, time_idx), 100);
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
    plot(r_x, r_y);
    % plot(t_seq, sqrt(v_x.^2 + v_y.^2));
    % poly_roots = roots(conv(v_x_poly, a_x_poly) + conv(v_y_poly, a_y_poly));
    % for j = 1 : length(poly_roots)
    %   imag_part = imag(poly_roots(j));
    %   real_part = real(poly_roots(j));
    %   if ~imag_part && real_part >= keyframe_list(i, 3) && real_part <= keyframe_list(i + 1, 3)
    %     real_roots = [real_roots, real_part];
    %   endif
    % endfor
    % low_cnt = 0;
    % high_cnt = 0;
    % for j = 1 : length(real_roots)
    %   vel_at_root = vel(real_roots(j));
    %   if vel_at_root < 0.5
    %     low_cnt += 1;
    %   endif
    %   if vel_at_root > 1
    %     high_cnt += 1;
    %   endif
    % endfor
  endfor
  % real_roots
  hold off
endfunction
