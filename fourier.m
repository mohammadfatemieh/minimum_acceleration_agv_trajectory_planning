classdef fourier
  % fourier  a class implement a basic Fourier-Euler basis of function on an interval
  %
  %   PROPERTIES
  %     order: order of the basis
  %     boundary: [left: double, right: double]
  %     series: orderx2 cell, contains trigonometric function or double
  %     amplitude: return 1x(1 + 2 * order) array, each denote the amplitude of a base
  %
  %   METHODS
  %     base(n, type): get the specified base with order n and type (0 -> cos, 1 -> sin)
  %     square_int_value
  properties (GetAccess=public, SetAcess=public)
    order
    boundary
    series
  endproperties
  properties (Dependent)
    amplitude
  endproperties
  methods
    function obj = fourier(order, boundary) % boundary = [left, right]
      obj.order = order;
      obj.boundary = boundary;
      obj.boundary(2) = 4 * boundary(2) - 3 * boundary(1);
      obj.series = cell(order + 1, 2);
      obj.series{1, 1} = sqrt(1 / diff(obj.boundary));
      obj.series{1, 2} = 0;
      amplitude = sqrt(2) * obj.series{1, 1};
      vertical_shift = 0;
      for n = 1 : order
        for j = 1 : 2
          angular_frequency = 2 * pi * n / diff(obj.boundary);
          phase_shift = -sum(obj.boundary) * pi * n / diff(obj.boundary);
          func_type_indicator = (j - 1) + (2 - j) * 1i;
          obj.series{n + 1, j} = trigonometric(func_type_indicator, amplitude, angular_frequency, phase_shift, vertical_shift);
        endfor
      endfor
    endfunction
    function amp = get.amplitude(obj)
      amp = zeros(1, 2 * obj.order + 1);
      amp(1) = 0;
      for n = 1 : obj.order
        for j = 1 : 2
          amp(n * 2 - 1 + j) = obj.series{n + 1, j}.amplitude;
        endfor
      endfor
    endfunction
    function b = base(obj, order, func_type_indicator)
      b = obj.series{order + 1, func_type_indicator + 1};
    endfunction
    function val = value(obj, x)
      val = zeros(1, 2 * obj.order + 1);
      val(1) = obj.series{1, 1};
      for n = 1 : obj.order
        for j = 1 : 2
          val(n * 2 - 1 + j) = obj.series{n + 1, j}.value(x);
        endfor
      endfor
    endfunction
    function val = square_int_value(obj, cur_x, nxt_x)
      % square_int_value  First calculate the Descartes product of two same basis,
      %                   then integrate it with know lower and upper limit
      %
      %   val = square_int_value(cur_x, nxt_x)
      %                         1 cos(x) sin(x) cos(2x) sin(2x)
      %            ------------------------------------------
      %     1      |            0 0      0      0       0      
      %     cos(x) |            0 c^2(x) 0      0       0      
      %     sin(x) | integrate( 0 0      s^2(x) 0       0      ) from cur_x to nxt_x
      %     cos(2x)|            0 0      0      c^2(2x) 0      
      %     sin(2x)|            0 0      0      0       s^2(2x)
      %     
      %     any only return the diagonal.
      val = zeros(1, 2 * obj.order + 1);
      val(1) = obj.series{1, 1};
      for n = 1 : obj.order
        for j = 1 : 2
          sign = 1;
          if obj.series{n + 1, j}.is_sin
            sign = -1;
          endif
          amplitude = obj.series{n + 1, j}.amplitude;
          angular_frequency = obj.series{n + 1, j}.angular_frequency;
          phase_shift = obj.series{n + 1, j}.phase_shift;
          vertical_shift = obj.series{n + 1, j}.vertical_shift;
          ans_trig = trigonometric(sign * (1 + 0i), ...
                                   amplitude^2 / (4 * angular_frequency), ...
                                   2 * angular_frequency, ...
                                   2 * phase_shift, ...
                                   amplitude^2 .* [0.5, phase_shift / (2 * angular_frequency)]);
          val(n * 2 - 1 + j) = ans_trig.value(nxt_x) - ans_trig.value(cur_x);
        endfor
      endfor
    endfunction
    function val = mul_der_int_value(obj, cur_t, nxt_t)
      der = obj.derivative;
      val = 0;
      for n = 1 : obj.order
        val_blk = zeros(2, 2);
        for j = 1 : 2
          ans_sign = obj.series{n + 1, j}.sign * der.series{n + 1, 3 - j}.sign;
          if obj.series{n + 1, j}.is_sin
            ans_sign *= -1;
          endif
          angular_frequency = 2 * obj.series{n + 1, j}.angular_frequency;
          phase_shift = 2 * obj.series{n + 1, j}.phase_shift;
          amplitude = obj.series{n + 1, j}.amplitude * der.series{n + 1, 3 - j}.amplitude / (4 * obj.series{n + 1, j}.angular_frequency);
          vertical_shift = ...
            obj.series{n + 1, j}.amplitude * der.series{n + 1, 3 - j}.amplitude * ...
            obj.series{n + 1, j}.sign * der.series{n + 1, 3 - j}.sign * ...
            [0.5, obj.series{n + 1, j}.phase_shift / (2 * obj.series{n +1, j}.angular_frequency)];
          ans_trig = trigonometric(ans_sign * (1 + 0i), ...
                                   amplitude, ...
                                   angular_frequency, ...
                                   phase_shift, ...
                                   vertical_shift);
          disp(obj.series{n + 1, j});
          disp(der.series{n + 1, 3 - j});
          disp(ans_trig);
          val_blk(j, 3 - j) = ans_trig.value(nxt_t) - ans_trig.value(cur_t);
        endfor
        val = blkdiag(val, val_blk);
      endfor
    endfunction
    function der_obj = derivative(obj)
      der_obj = obj;
      der_obj.series{1, 1} = 0;
      der_obj.series{1, 2} = 0;
      for n = 1 : obj.order
        for j = 1 : 2
          der_obj.series{n + 1, j} = obj.series{n + 1, j}.derivative;
        endfor
      endfor
    endfunction
    function scale_obj = scale(obj, coeff)
      % scale  Change the amplitude of each one base.
      %
      %   scaled_fourier = f.scale(coeff)
      %     coef  an array with size 1x(1 + 2 * order)
      %     Return a fourier object with each base have its amplitude multiplied
      %     by a number in coeff. This number is specified as
      %       [1, cos(x), sin(x), cos(2x), sin(2x), ...]
      scale_obj = obj;
      scale_obj.series{1, 1} *= coeff(1);
      for n = 1 : obj.order
        for j = 1 : 2
          scale_obj.series{n + 1, j} = scale_obj.series{n + 1, j}.change_amplitude(scale_obj.series{n + 1, j}.amplitude * coeff(n * 2 - 1 + j));
        endfor
      endfor
    endfunction
    function disp(obj)
      disp(obj.series{1, 1})
      for n = 1 : obj.order
        for j = 1 : 2
          disp(obj.series{n + 1, j});
        endfor
      endfor
    endfunction
    function plot(obj)
      hold on
      for n = 1 : obj.order
        for j = 1 : 2
          plot(obj.series{n + 1, j});
        endfor
      endfor
      hold off
    endfunction
  endmethods
endclassdef
