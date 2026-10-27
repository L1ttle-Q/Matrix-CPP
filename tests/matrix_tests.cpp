#define MATRIX_CPP_NO_MAIN
#include "../matrix.cpp"

int main()
{
    int failures = 0;
    auto check = [&](bool success, const char* name) {
        std::cout << (success ? "PASS " : "FAIL ") << name << '\n';
        if (!success) ++failures;
    };
    Matrix<double> a(2, 3), b(2, 3);
    b[2][3] = -0.75;
    check(std::abs(norm(a, b) - 0.75) < 1e-12, "rectangular norm includes last column and fractions");
    Matrix<double> x(3, 1), y(3, 1);
    y[3][1] = 0.125;
    check(std::abs(norm(x, y) - 0.125) < 1e-12, "column vector norm");
    check(norm(a, x) == -1, "norm rejects incompatible shapes");
    Matrix<int> m(2, 2);
    m[1][1] = 4; m[1][2] = 7; m[2][1] = 2; m[2][2] = 6;
    check(m * Matrix<int>::I(2) == m, "identity multiplication");
    Matrix<int> transposed = m.T();
    check(transposed[1][2] == 2 && transposed[2][1] == 7, "transpose");
    check(norm(m * m.inverse(), Matrix<double>::I(2)) < 1e-10, "inverse residual");
    Matrix<int> copy(m);
    copy[1][1] = -5;
    check(m[1][1] == 4, "copy owns independent storage");
    return failures ? 1 : 0;
}
