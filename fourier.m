classdef fourier
  properties (GetAccess=public, SetAcess=private)
    order
    boundary
    series
  endproperties
  methods
    function obj = fourier(order, boundary) % boundary = [left, right]
      obj.order = order;
      obj.boundary = boundary;
      obj.series = cell(order + 1, 2);
      obj.series{1, 1} = sqrt(1 / diff(boundary));
      obj.series{1, 2} = 0;
      amplitude = sqrt(2) * obj.series{1, 1};
      vertical_shift = 0;
      for n = 1 : order
        for j = 1 : 2
          angular_frequency = 2 * pi * n / diff(boundary);
          phase_shift = -sum(boundary) * pi * n / diff(boundary);
          obj.series{n + 1, j} = trigonometric(j - 1, amplitude, angular_frequency, phase_shift, vertical_shift);
        endfor
      endfor
    endfunction
    function b = base(obj, order, is_sin)
      b = obj.series{order + 1, is_sin + 1};
    endfunction
    function val = value(obj, x)
      val = zeros(1, 2 * order + 1);
      val(1) = obj.series{1, 1};
      for n = 1 : order
        for j = 1 : 2
          val(n * 2 - 1 + j) = obj.series{n + 1, j}.value(x);
        endfor
      endfor
    endfunction
    function val = square_int_value(obj, cur_x, nxt_x)
      val = zeros(1, 2 * obj.order + 1);
      val(1) = obj.series{1, 1};
      for n = 1 : obj.order
        for j = 1 : 2
          obj.series{n + 1, j}
          amplitude = obj.series{n + 1, j}.amplitude;
          angular_frequency = obj.series{n + 1, j}.angular_frequency;
          phase_shift = obj.series{n + 1, j}.phase_shift;
          vertical_shift = obj.series{n + 1, j}.vertical_shift;
          ans_trig = trigonometric(1, ...
                                   (-1)^(j) * amplitude^2 / (4 * angular_frequency), ...
                                   - 2 * angular_frequency, ...
                                   2 * phase_shift, ...
                                   amplitude^2 .* [0.5, -phase_shift / (2 * angular_frequency)])
          val(n * 2 - 1 + j) = ans_trig.value(nxt_x) - ans_trig.value(cur_x)
        endfor
      endfor
    endfunction
    function derivative(obj)
    endfunction
    function integral(obj)
    endfunction
    function mut_der(obj)
    endfunction
    function mut_int(obj)
    endfunction
  endmethods
endclassdef
