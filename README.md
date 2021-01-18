# Minimum acceleration trajectory planning for AGV

In reference[1], the author model this problem as a quadratic programming problem,
with integral of square of snap over time as the target function.

In this case, we use the Fourier-Euler as the basis, rather than polynomials,
to improve the instability of the original solution when the order of the
polynomial gets higher.

## Mathematical Notion

Suppose we have `m` keyframes and thus `m - 1` segments. We need one function for
each segment. Denote `pos_m` as the function of position for the `m`th segment,
we could write

    pos_m(t) = a_1_m + a_2_m * cos(t) + a_3_m * sin(t) + a_4_m * cos(2t) + a_5_m * cos(2t) ...
             = [1, cos(t), sin(t), cos(2t), sin(2t), ...] * [a_1_m, a_2_m, a_3_m, a_4_m, a_5_m, ...]'

If we use first 5 base of the Fourier-Euler basis, we could have each segment
represented as the product of a 1x5 and a 5x1 matrix. In total we have `5 * m`
coefficients, which make up `x`.

In this project, we denote the cost function as

    0.5 * x * H * x + f * x

And the constrain as

    A * x = b

## Reference

[1] D. Mellinger, and V. Kumar, "Minimum Snap Trajectory Generation and Control for
Quadrotors," *IEEE International Conference on Robotics and Automation*,
Shanghai, China, May 9-13, 2011.
