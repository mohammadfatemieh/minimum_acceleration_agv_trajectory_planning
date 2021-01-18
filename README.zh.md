# AGV 最小加速度轨迹生成

本项目基于论文[1] 实现。在该论文的基础上，本项目使用傅里叶-欧拉级数代替了论文中的多项式作为基函数，来改善原文方法在高阶多项式拟合的情况下的不稳定性。由于 AGV 和四轴飞行器的运动学差异，本项目的优化目标为最小加速度。

## 依赖

带有 optim 包的 GNU Octave 或者 MATLAB 。

## 使用

### GNU Octave

在 GNU Octave 的终端交互界面加载 optim 包，然后调用 `ui.m` ，并输入相应参数

```octave
octave:1> pkg load optim
octave:2> ui
geometry (e.g. [x, y]): % 这里输入机器人的尺寸
maximum velocity: % 这里输入机器人的最大速度
theta at the start point (in degree): % 这里以角度输入起点机器人的位姿
theta at the end point (in degree): % 这里以角度输入终点机器人的位姿
velocity at the start point (magnitude): % 这里输入起点机器人的速度大小
velocity at the end point (magnitude): % 这里输入终点机器人的速度大小
```

按下回车后，会有新窗口弹出，使用鼠标在弹出窗口的座标轴内点击一下设置障碍物的一端，再点击一次设置障碍物的另一端。点击设置机器人的起点，再次点击设置机器人的终点。等待至窗口内有轨迹生成，然后在 GNU Octave 的终端交互环境中按下回车，绘图窗口中会以蓝色圆圈指示当前小车的质心位置，多次按回车可追踪小车的运动。在小车的运动过程中，以鼠标左键点击窗口，可设置新的终点并触发重新规划。当小车完成一段时间的运动，或重新设置终点后，会有窗口绘出小车的速度图像。

### MATLAB（未测试）

进入本项目的文件夹。在 MATLAB 下方的命令窗口输入 `ui` ，然后参见 GNU Octave 一节。

## 一些想法

以傅里叶级数来拟合机器人的运动轨迹，一是因为机器人的运动速度有上限，二是考虑到快速傅里叶变换可能带来的性能提升。虽然最后并没有用上快速傅里叶变换，但是如果不采使用 MATLAB 提供的二次规划求解器，而是自己实现，那么在二次规划的过程中使用快速傅里叶变化可能可以减少矩阵的规模，从而降低内存消耗，提高性能。

## 参考资料

[1] D. Mellinger, and V. Kumar, "Minimum Snap Trajectory Generation and Control for
Quadrotors," *IEEE International Conference on Robotics and Automation*,
Shanghai, China, May 9-13, 2011.
[2] [Introduction to Object-Oriented Programming in MATLAB](https://www.mathworks.com/company/newsletters/articles/introduction-to-object-oriented-programming-in-matlab.html)
[3] [Limited-memory BFGS - Wikipedia](https://en.wikipedia.org/wiki/Limited-memory_BFGS)
