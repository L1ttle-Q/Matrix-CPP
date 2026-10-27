# Matrix-CPP

单文件 C++11 矩阵模板与数值迭代学习项目。`matrix.cpp` 同时包含 `Matrix<T>` 实现和一个固定 5×5 线性方程组演示。

## 构建与运行

需要 Git、CMake 3.16+ 和支持 C++11 的编译器，无第三方库或 Python 依赖。

Windows + MinGW-w64：

```powershell
git clone https://github.com/L1ttle-Q/Matrix-CPP.git
cd Matrix-CPP
cmake -S . -B build -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
.\build\matrix_demo.exe
ctest --test-dir build --output-on-failure
```

Linux/macOS 配置省略 `-G "MinGW Makefiles"`，运行 `./build/matrix_demo`。Visual Studio 使用单独的 `build-vs` 目录，构建加 `--config Release`，运行 `build-vs/Release/matrix_demo.exe`，CTest 加 `-C Release`。

也可直接编译演示：

```sh
g++ -std=c++11 -O2 matrix.cpp -o matrix_demo
```

Windows 对应输出可命名 `matrix_demo.exe`。程序无需命令行参数或标准输入；矩阵和右端向量在 `main()` 中定义。输出依次为 Jacobi、Gauss-Seidel、SOR 的迭代过程和逆矩阵求解结果。默认最多迭代 100 轮，阈值为 `1e-6`；演示矩阵不保证迭代收敛，`non-convergent!` 不是编译或启动失败。

## 在自己的代码中使用

当前模板未拆成独立头文件；可用宏屏蔽演示入口后包含一次实现文件：

```cpp
#define MATRIX_CPP_NO_MAIN
#include "matrix.cpp"

int main()
{
    Matrix<double> a(2, 2);
    a[1][1] = 4; a[1][2] = 7;
    a[2][1] = 2; a[2][2] = 6;
    (a * a.inverse()).Print(); // 近似 2×2 单位矩阵
}
```

保存为仓库根目录下的 `example.cpp`，编译 `g++ -std=c++11 example.cpp -o example`；不要同时把 `matrix.cpp` 作为第二个编译单元传入。当前这种包含方式适用于单编译单元示例。

## 接口速查

| 接口 | 说明 |
| --- | --- |
| `Matrix<T>(rows, cols)` | 创建零初始化矩阵；默认 100×100 |
| `a[i][j]` | **从 1 开始**的行列下标，不检查越界 |
| `get_r()` / `get_c()` | 行数 / 列数 |
| `+`、`-`、`*`、`==` | 矩阵运算，另支持右乘 `double` 标量 |
| `Matrix<T>::I(n)` / `T()` | 单位矩阵 / 转置 |
| `inverse()` | Gauss-Jordan 逆矩阵，返回 `Matrix<double>` |
| `subMatrix(x1,y1,x2,y2)` | 按闭区间取子矩阵；现有边界处理有局限 |
| `norm(a,b)` | 最大元素绝对差；形状不一致返回 `-1` |
| `Jacobi` / `Gauss_Seidel` / `SOR` | 单步迭代函数；`Solve` 执行并打印迭代 |

## 测试与数值限制

`tests/matrix_tests.cpp` 有 7 项检查：非方阵与列向量误差、形状不匹配、单位矩阵乘法、转置、求逆残差及深复制。它们注册为一项 CTest，失败时返回非零，Release 构建也不会移除检查。

本次修复了 `norm` 把行数当列数的越界问题，并显式用 `double` 的绝对值比较，保留小数误差。

`inverse()` 用精确零判断主元，没有数值容差；`subMatrix()` 不支持所有单行/单列范围，且部分边界缺少检查；部分失败运算返回 `NullMatrix` 而不抛异常。SOR 的松弛系数当前固定为 `1.22`。测试未覆盖全部数值稳定性和无效维度，项目适合学习与实验。

2026-10-26：Windows、MinGW-w64 GCC 8.1、CMake 3.31 下，演示构建及含 7 项检查的 CTest 通过。其他系统未在本轮实测。

使用 [Apache License 2.0](LICENSE)。
