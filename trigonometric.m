classdef trigonometric
  properties (GetAccess=public, SetAcess=private)
    func_type_indicator
    amplitude
    angular_frequency
    phase_shift
    vertical_shift
  endproperties
  properties (Dependent)
    period
    sign
    func
    is_sin
    is_cos
  endproperties
  methods
    function f = trigonometric(func_type_indicator, amplitude, angular_frequency, phase_shift, vertical_shift)
      f.func_type_indicator = func_type_indicator;
      f.amplitude = amplitude;
      f.angular_frequency = angular_frequency;
      f.phase_shift = phase_shift;
      f.vertical_shift = vertical_shift;
    endfunction
    function period = get.period(obj)
      period = 2 * pi / obj.angular_frequency;
    endfunction
    function f = change_amplitude(obj, new_amplitude)
      f = trigonometric(obj.func_type_indicator, ...
                        new_amplitude, ...
                        obj.angular_frequency, ...
                        obj.phase_shift, ...
                        obj.vertical_shift);
    endfunction
    function sign = get.sign(obj)
      sign = 1;
      if obj.func_type_indicator == -1 + 0i || obj.func_type_indicator == 0 - 1i
        sign = -1;
      endif
    endfunction
    function func = get.func(obj)
      func = @(x) sin(x);
      if obj.func_type_indicator == 0 + 1i || obj.func_type_indicator == 0 - 1i
        func = @(x) cos(x);
      endif
    endfunction
    function ret = get.is_sin(obj)
      if obj.func_type_indicator == 1 + 0i || obj.func_type_indicator == -1 + 0i
        ret = logical(1);
      else
        ret = logical(0);
      endif
    endfunction
    function ret = get.is_cos(obj)
      ret = ~obj.is_sin;
    endfunction
    function val = value(obj, x)
      val = obj.sign * obj.amplitude * obj.func(obj.angular_frequency * x + obj.phase_shift) + polyval(obj.vertical_shift, x);
    endfunction
    function der_trig = derivative(obj)
      der_trig = trigonometric(obj.func_type_indicator * 1i, ...
                               obj.amplitude * obj.angular_frequency, ...
                               obj.angular_frequency, ...
                               obj.phase_shift, ...
                               polyder(obj.vertical_shift));
    endfunction
    function int_trig = integral(obj)
      int_trig = trigonometric(obj.func_type_indicator / 1i, ...
                               obj.amplitude / obj.angular_frequency, ...
                               obj.angular_frequency, ...
                               obj.phase_shift, ...
                               polyint(obj.vertical_shift));
    endfunction
    function plot(obj)
      x = linspace(0, 1.75 * obj.period, 100);
      y = obj.value(x);
      plot(x, y);
    endfunction
    function disp(obj)
      if obj.amplitude == 0
        output_str = '0';
      else
        if obj.func_type_indicator == 1 + 0i
          output_str = [num2str(obj.amplitude), ' sin(', num2str(obj.angular_frequency), ' x'];
        else
          if obj.func_type_indicator == 0 + 1i
            output_str = [num2str(obj.amplitude), ' cos(', num2str(obj.angular_frequency), ' x'];
          else
            if obj.func_type_indicator == -1 + 0i
              output_str = [num2str(-obj.amplitude), ' sin(', num2str(obj.angular_frequency), ' x'];
            else % obj.func_type_indicator == 0 - 1i
              output_str = [num2str(-obj.amplitude), ' cos(', num2str(obj.angular_frequency), ' x'];
            endif
          endif
        endif
        if obj.phase_shift ~= 0
          output_str = [output_str, ' + ', num2str(obj.phase_shift), ')'];
        else
          output_str = [output_str, ')'];
        endif
        if obj.vertical_shift ~= 0
          output_str = [output_str, ' + ', num2str(obj.vertical_shift)];
        endif
      endif
      disp(output_str);
    endfunction
  endmethods
endclassdef
