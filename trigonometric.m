classdef trigonometric
  properties (GetAccess=public, SetAcess=private)
    is_sin
    amplitude
    angular_frequency
    phase_shift
    vertical_shift
  endproperties
  properties (Dependent)
    period
  endproperties
  methods
    function f = trigonometric(is_sin, amplitude, angular_frequency, phase_shift, vertical_shift)
      f.is_sin = is_sin;
      f.amplitude = amplitude;
      f.angular_frequency = angular_frequency;
      f.phase_shift = phase_shift;
      f.vertical_shift = vertical_shift;
    endfunction
    function period = get.period(obj)
      period = 2 * pi / obj.angular_frequency;
    endfunction
    function f = change_amplitude(obj, new_amplitude)
      f = trigonometric(obj.is_sin, ...
                        new_amplitude, ...
                        obj.angular_frequency, ...
                        obj.phase_shift, ...
                        obj.vertical_shift);
    endfunction
    function val = value(obj, x)
      if obj.is_sin == 1 + 0i
        val = obj.amplitude * sin(obj.angular_frequency * x + obj.phase_shift) + polyval(obj.vertical_shift, x);
      else
        if obj.is_sin == 0 + 1i
          val = obj.amplitude * cos(obj.angular_frequency * x + obj.phase_shift) + polyval(obj.vertical_shift, x);
        else
          if obj.is_sin == -1 + 0i
            val = -obj.amplitude * sin(obj.angular_frequency * x + obj.phase_shift) + polyval(obj.vertical_shift, x);
          else % obj.is_sin == 0 - 1i
            val = -obj.amplitude * cos(obj.angular_frequency * x + obj.phase_shift) + polyval(obj.vertical_shift, x);
          endif
        endif
      endif
    endfunction
    function der_trig = derivative(obj)
      if obj.is_sin
        der_trig = trigonometric(obj.is_sin * i, ...
                                 obj.amplitude * obj.angular_frequency, ...
                                 obj.angular_frequency, ...
                                 obj.phase_shift, ...
                                 polyder(obj.vertical_shift));
      else
        der_trig = trigonometric(obj.is_sin * i, ...
                                 obj.amplitude * obj.angular_frequency, ...
                                 obj.angular_frequency, ...
                                 obj.phase_shift, ...
                                 polyder(obj.vertical_shift));
      endif
    endfunction
    function int_trig = integral(obj)
      if obj.is_sin
        int_trig = trigonometric(~obj.is_sin, ...
                                 -obj.amplitude / obj.angular_frequency, ...
                                 obj.angular_frequency, ...
                                 obj.phase_shift, ...
                                 polyint(obj.vertical_shift));
      else
        int_trig = trigonometric(~obj.is_sin, ...
                                 obj.amplitude / obj.angular_frequency, ...
                                 obj.angular_frequency, ...
                                 obj.phase_shift, ...
                                 polyint(obj.vertical_shift));
      endif
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
        if obj.amplitude ~= 1
          output_str = [num2str(obj.amplitude), ' '];
        else
          output_str = '';
        endif
        if obj.is_sin
          output_str = [output_str, 'sin(', num2str(obj.angular_frequency), ' x'];
        else
          output_str = [output_str, 'cos(', num2str(obj.angular_frequency), ' x'];
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
