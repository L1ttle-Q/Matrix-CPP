#include <cstdio>
#include <type_traits>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <vector>
#include <cmath>

template<typename T1>
class Matrix
{
private:
    int r, c;
    T1 **matrix;
public:
    static const Matrix<T1> NullMatrix;

    Matrix(int n = 100, int m = 100): r(n), c(m)
    {
        if (!r || !c) {matrix = nullptr; return ;}
        matrix = new T1* [r + 1];
        for (int i = 1; i <= r; i++)
        {
            matrix[i] = new T1 [c + 1];
            for (int j = 1; j <= c; j++)
                matrix[i][j] = 0;
        }
    }

    template<typename U>
    Matrix(const Matrix<U> &b): r(b.get_r()), c(b.get_c())
    {
        if (!r || !c) {matrix = nullptr; return ;}
        matrix = new T1* [r + 1];
        for (int i = 1; i <= r; i++)
        {
            matrix[i] = new T1 [c + 1];
            for (int j = 1; j <= c; j++)
                matrix[i][j] = T1(b[i][j]);
        }
    }

    Matrix(const Matrix<T1> &b): r(b.get_r()), c(b.get_c())
    {
        if (!r || !c) {matrix = nullptr; return ;}
        matrix = new T1* [r + 1];
        for (int i = 1; i <= r; i++)
        {
            matrix[i] = new T1 [c + 1];
            for (int j = 1; j <= c; j++)
                matrix[i][j] = b[i][j];
        }
    }

    Matrix<double> inverse() const
    {
        if (!r || !c || r != c) return Matrix<double>::NullMatrix;

        int n = r;
        Matrix<double> aug(n, 2 * n);

        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= n; j++)
                aug[i][j] = matrix[i][j];

        for (int i = 1; i <= n; i++) aug[i][i + n] = 1;

        for (int i = 1; i <= n; i++)
        {
            int pivot = i;
            for (int j = i + 1; j <= n; j++)
                if (abs(aug[j][i]) > abs(aug[pivot][i]))
                    pivot = j;

            if (!aug[pivot][i]) return Matrix<double>::NullMatrix;

            if (i != pivot)
                for (int j = 1; j <= 2 * n; j++)
                    std::swap(aug[i][j], aug[pivot][j]);

            double pivot_val = aug[i][i];
            for (int j = 1; j <= 2 * n; j++)
                aug[i][j] /= pivot_val;

            for (int j = 1; j <= n; j++)
                if (j != i)
                {
                    double factor = aug[j][i];
                    for (int k = 1; k <= 2 * n; k++)
                        aug[j][k] -= factor * aug[i][k];
                }
        }

        return aug.subMatrix(1, n + 1, n, 2 * n);
    }

    Matrix T()
    {
        if (!r || !c) return NullMatrix;
        Matrix res(c, r);
        for (int i = 1; i <= r; i++)
            for (int j = 1; j <= c; j++)
                res[j][i] = matrix[i][j];
        return res;
    }

    ~Matrix()
    {
        for (int i = 1; i <= r;i++) delete[] matrix[i];
        delete[] matrix;
    }

    Matrix subMatrix(const int& x1, const int& y1, const int& x2, const int& y2) const
    {
        if (x2 < 0 || y2 < 0 || x1 >= r || y1 >= c) return NullMatrix;
        if (x2 <= x1 || y2 <= y1) return NullMatrix;
        Matrix res(x2 - x1 + 1, y2 - y1 + 1);
        for (int i = x1; i <= x2; i++)
            for (int j = y1; j <= y2; j++)
                res[i - x1 + 1][j - y1 + 1] = matrix[i][j];
        return res;
    }

    static Matrix I(int n)
    {
        Matrix I(n, n);
        for (int i = 1; i <= n; i++)
            I[i][i] = 1;
        return I;
    }

    Matrix& operator =(const Matrix<T1>& b)
    {
        if (*this == b) return *this;

        for (int i = 1; i <= r; i++)
            delete[] matrix[i];
        delete[] matrix;

        r = b.get_r();
        c = b.get_c();
        matrix = new T1* [r + 1];
        for (int i = 1; i <= r; i++)
        {
            matrix[i] = new T1 [c + 1];
            for (int j = 1; j <= c; j++)
                matrix[i][j] = b[i][j];
        }
        return *this;
    }

    template<typename U>
    Matrix<typename std::common_type<T1, U>::type> operator *(const Matrix<U> &b) const
    {
        using ResultType = typename std::common_type<T1, U>::type;
        if (c != b.get_r()) return Matrix<ResultType>::NullMatrix;
        Matrix<ResultType> res(r, b.get_c());
        for (int i = 1; i <= r; i++)
        {
            ResultType tmp;
            for (int k = 1; k <= c; k++)
            {
                tmp = matrix[i][k];
                for (int j = 1; j <= b.get_c(); j++)
                    res[i][j] += tmp * b[k][j];
            }
        }
        return res;
    }

    Matrix<T1> operator *(const Matrix<T1> &b) const
    {
        if (c != b.get_r()) return Matrix<T1>::NullMatrix;
        Matrix<T1> res(r, b.get_c());
        for (int i = 1; i <= r; i++)
        {
            T1 tmp;
            for (int k = 1; k <= c; k++)
            {
                tmp = matrix[i][k];
                for (int j = 1; j <= b.get_c(); j++)
                    res[i][j] += tmp * b[k][j];
            }
        }
        return res;
    }

    Matrix<double> operator *(const double &b) const
    {
        Matrix<double> res(r, c);
        for (int i = 1; i <= r; i++)
            for (int j = 1; j <= c; j++)
                res[i][j] = matrix[i][j] * b;
        return res;
    }

    template<typename U>
    Matrix<typename std::common_type<T1, U>::type> operator +(const Matrix<U> &b) const
    {
        using ResultType = typename std::common_type<T1, U>::type;
        if (r != b.get_r() || c != b.get_c()) return Matrix<ResultType>::NullMatrix;
        Matrix<ResultType> res(r, b.get_c());
        for (int i = 1; i <= r; i++)
            for (int j = 1; j <= c; j++)
                res[i][j] = ResultType(matrix[i][j]) + ResultType(b[i][j]);
        return res;
    }

    Matrix operator +(const Matrix& b) const
    {
        if (r != b.get_r() || c != b.get_c()) return Matrix<T1>::NullMatrix;
        Matrix res(r, c);
        for (int i = 1; i <= r; i++)
            for (int j = 1; j <= c; j++)
                res[i][j] = matrix[i][j] + b[i][j];
        return res;
    }

    template<typename U>
    Matrix<typename std::common_type<T1, U>::type> operator -(const Matrix<U> &b) const
    {
        using ResultType = typename std::common_type<T1, U>::type;
        if (r != b.get_r() || c != b.get_c()) return Matrix<ResultType>::NullMatrix;
        Matrix<ResultType> res(r, b.get_c());
        for (int i = 1; i <= r; i++)
            for (int j = 1; j <= c; j++)
                res[i][j] = ResultType(matrix[i][j]) - ResultType(b[i][j]);
        return res;
    }

    Matrix operator -(const Matrix& b) const
    {
        if (r != b.get_r() || c != b.get_c()) return Matrix<T1>::NullMatrix;
        Matrix res(r, c);
        for (int i = 1; i <= r; i++)
            for (int j = 1; j <= c; j++)
                res[i][j] = matrix[i][j] - b[i][j];
        return res;
    }

    bool operator ==(const Matrix& b) const
    {
        if (r != b.get_r() || c != b.get_c()) return false;
        for (int i = 1; i <= r; i++)
            for (int j = 1; j <= c; j++)
                if (matrix[i][j] != b[i][j])
                    return false;
        return true;
    }

    T1* operator [](int x) const
    {
        return matrix[x];
    }

    void Print() const
    {
        for (int i = 1; i <= r; i++)
        {
            if (i == 1) std::cout << "[";
            std::cout << "[";
            for (int j = 1; j < c; j++)
                std::cout << matrix[i][j] << " ";
            std::cout << matrix[i][c] << "]";
            if (i == r) std::cout << "]";
            std::cout << std::endl;
        }
    }

    int get_r() const {return r;}
    int get_c() const {return c;}
};

template<>
void Matrix<double>::Print() const
{
    for (int i = 1; i <= r; i++)
    {
        if (i == 1) putchar('[');
        putchar('[');
        for (int j = 1; j < c; j++)
            printf("%.8f ", matrix[i][j]);
        printf("%.8f]", matrix[i][c]);
        if (i == r) putchar(']');
        putchar('\n');
    }
}

template<>
void Matrix<int>::Print() const
{
    for (int i = 1; i <= r; i++)
    {
        if (i == 1) putchar('[');
        putchar('[');
        for (int j = 1; j < c; j++)
            printf("%d ", matrix[i][j]);
        printf("%d]", matrix[i][c]);
        if (i == r) putchar(']');
        putchar('\n');
    }
}

template<typename T1, typename T2>
double norm(const Matrix<T1>& a, const Matrix<T2>& b)
{
    if (a.get_r() != b.get_r() || a.get_c() != b.get_c()) return -1;

    // double res = 0;
    // for (int i = 1; i <= a.get_r(); i++)
    //     for (int j = 1; j <= a.get_c(); j++)
    //         res += pow(a[i][j] - b[i][j], 2);
    // return res / (a.get_r() * a.get_c());
    double res = 0;
    for (int i = 1; i <= a.get_r(); i++)
        for (int j = 1; j <= a.get_r(); j++)
            res = std::max(res, abs(a[i][j] - b[i][j]));
    return res;
}

template<typename T1>
const Matrix<T1> Matrix<T1>::NullMatrix = Matrix(0, 0);

template<typename T>
Matrix<T> get_U(Matrix<T> mat)
{
    if (!mat.get_c() || mat.get_c() != mat.get_r())
        return Matrix<T>::NullMatrix;

    int n = mat.get_c();
    Matrix<T> U(mat);
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= i; j++)
            U[i][j] = 0;
    return U;
}

template<typename T>
Matrix<T> get_D(Matrix<T> mat)
{
    if (!mat.get_c() || mat.get_c() != mat.get_r())
        return Matrix<T>::NullMatrix;

    int n = mat.get_c();
    Matrix<T> D(mat);
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            if (i != j)
                D[i][j] = 0;
    return D;
}

template<typename T>
Matrix<T> get_L(Matrix<T> mat)
{
    if (!mat.get_c() || mat.get_c() != mat.get_r())
        return Matrix<T>::NullMatrix;

    int n = mat.get_c();
    Matrix<T> L(mat);
    for (int i = 1; i <= n; i++)
        for (int j = i; j <= n; j++)
            L[i][j] = 0;
    return L;
}

template<typename T1, typename T2>
Matrix<double> Jacobi(Matrix<double> x, Matrix<T1> A, Matrix<T2> b)
{
    Matrix<double> U(get_U(A)), L(get_L(A)), D(get_D(A));
    Matrix<double> B((D.inverse() * (L + U)) * (-1.0)), g(D.inverse() * b);
    return B * x + g;
}

template<typename T1, typename T2>
Matrix<double> Gauss_Seidel(Matrix<double> x, Matrix<T1> A, Matrix<T2> b)
{
    Matrix<double> U(get_U(A)), L(get_L(A)), D(get_D(A));
    Matrix<double> B((D + L).inverse() * (-1.0) * U), g((D + L).inverse() * b);
    return B * x + g;
}

template<typename T1, typename T2>
Matrix<double> SOR(Matrix<double> x, Matrix<T1> A, Matrix<T2> b)
{
    double w = 1.22; // 1.8 1.22
    // Matrix<double> I = Matrix<double>::I(A.get_r());
    // Matrix<double> U(get_U(A)), L(get_L(A)), D(get_D(A));
    // Matrix<double> B((I + D.inverse() * L * w).inverse() * (I * (1 - w) - D.inverse() * U * w)),
    //                g((I + D.inverse() * L * w) * D.inverse() * b * w);
    // Matrix<double> B((D + L).inverse() * (-1.0) * U), g((D + L).inverse() * b);
    return Gauss_Seidel(x, A, b) * w + x * (1 - w);
}

template<typename T1, typename T2, typename T3>
void Solve(Matrix<T1> init_x, Matrix<T2> A, Matrix<T3> b,
           Matrix<double> (*kernal)(Matrix<double> x, Matrix<T2> A, Matrix<T3> b),
           int epoch = 100, double eps = 1e-6)
{
    Matrix<double> x(init_x);
    for (int i = 1; i <= epoch; i++)
    {
        printf("epoch %d:\n", i);
        Matrix<double> tmp(kernal(x, A, b));
        tmp.Print();
        if (norm(x, tmp) < eps){printf("convergent!\n");return ;}
        x = tmp;
    }
    printf("non-convergent!\n");
}

int main()
{
    // Matrix<double> init_x(3, 1);
    // Matrix<double> A(3, 3);
    // Matrix<double> b(3, 1);
    // init_x[1][1] = 1; init_x[2][1] = 1; init_x[3][1] = 1;

    // A[1][1] = 4; A[1][2] = 3; A[1][3] = 0;
    // A[2][1] = 3; A[2][2] = 4; A[2][3] = -1;
    // A[3][1] = 0; A[3][2] = -1; A[3][3] = 4;
    // b[1][1] = 24; b[2][1] = 30; b[3][1] = -24;

    // printf("Jacobi:\n");
    // Solve(init_x, A, b, Jacobi);
    // printf("\nGauss-Seidel:\n");
    // Solve(init_x, A, b, Gauss_Seidel);
    // printf("\nSOR:\n");
    // Solve(init_x, A, b, SOR);

    // printf("\nGaussian Elimination:\n");
    // (A.inverse() * b).Print();

    Matrix<double> init_x(5, 1);
    Matrix<double> A(5, 5);
    Matrix<double> b(5, 1);
    
    for (int i = 1; i <= 5; ++i) {
        init_x[i][1] = 1;
    }

    A[1][1] = 0.8147; A[1][2] = 0.0975; A[1][3] = 0.1576; A[1][4] = 0.1419; A[1][5] = 0.6557;
    A[2][1] = 0.9058; A[2][2] = 0.2785; A[2][3] = 0.9706; A[2][4] = 0.4218; A[2][5] = 0.0357;
    A[3][1] = 0.1270 * 1e10; A[3][2] = 0.5469; A[3][3] = 0.9572; A[3][4] = 0.9157; A[3][5] = 0.8491;
    A[4][1] = 0.9134; A[4][2] = 0.9575; A[4][3] = 0.4854 * 1e8; A[4][4] = 0.7922; A[4][5] = 0.9340;
    A[5][1] = 0.6324; A[5][2] = 0.9649; A[5][3] = 0.8003; A[5][4] = 0.9595; A[5][5] = 0.6787;

    b[1][1] = 2.258000e-9 * 1e9;
    b[2][1] = 1.597700e-9 * 1e9;
    b[3][1] = 1.270000002354900 * 1e9;
    b[4][1] = 0.024270003904200 * 1e9;
    b[5][1] = 3.360250e-9 * 1e9;

    printf("Jacobi:\n");
    Solve(init_x, A, b, Jacobi);
    printf("\nGauss-Seidel:\n");
    Solve(init_x, A, b, Gauss_Seidel);
    printf("\nSOR:\n");
    Solve(init_x, A, b, SOR);

    printf("\nGaussian Elimination:\n");
    (A.inverse() * b).Print();

    return 0;
}