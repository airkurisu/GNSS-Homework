/*
 * src/matrix.cpp - 矩阵运算具体实现 (Row-Major Version)
 * 项目: MYGNSS10
 */

#include "matrix.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>

using namespace std;

// 矩阵乘法: C = alpha * op(A) * op(B) + beta * C
// 实现说明: 采用行主序 (Row-Major) 索引: index = row * num_cols + col
void matmul(const char *tr, int n, int k, int m, double alpha,
            const double *A, const double *B, double beta, double *C)
{
    double d;
    // 解析转置标志: f=1(NN), f=2(NT), f=3(TN), f=4(TT)
    int f = tr[0] == 'N' ? (tr[1] == 'N' ? 1 : 2) : (tr[1] == 'N' ? 3 : 4);

    for (int i = 0; i < n; i++) {       // i: C 的行索引
        for (int j = 0; j < m; j++) {   // j: C 的列索引
            d = 0.0;
            for (int x = 0; x < k; x++) { // x: 公共维度索引

                // --- 获取 A 的元素 op(A)[i][x] ---
                // 若 A 不转置 (N): A 是 n x k 矩阵, 元素为 A[i][x] -> A[i*k + x]
                // 若 A 转置 (T)  : A 是 k x n 矩阵, 元素为 A[x][i] -> A[x*n + i]
                double a = (f == 1 || f == 2) ? A[i * k + x] : A[x * n + i];

                // --- 获取 B 的元素 op(B)[x][j] ---
                // 若 B 不转置 (N): B 是 k x m 矩阵, 元素为 B[x][j] -> B[x*m + j]
                // 若 B 转置 (T)  : B 是 m x k 矩阵, 元素为 B[j][x] -> B[j*k + x]
                double b = (f == 1 || f == 3) ? B[x * m + j] : B[j * k + x];

                d += a * b;
            }

            // 计算结果 C[i][j] -> C[i*m + j]
            if (beta == 0.0) C[i * m + j] = alpha * d;
            else C[i * m + j] = alpha * d + beta * C[i * m + j];
        }
    }
}

// 矩阵求逆 (Gauss-Jordan, Row-Major)
int matinv(double *A, int n) {
    if (n < 1) return -1;
    double scale = 0;
    for (int q=0; q<n*n; ++q) {
        if (!std::isfinite(A[q])) return -1;
        scale = std::max(scale, std::fabs(A[q]));
    }
    if (scale == 0) return -1;
    int i, j, k;
    // 记录行交换、列交换和主元状态
    vector<int> indxc(n), indxr(n), ipiv(n, 0);

    for (i = 0; i < n; i++) {
        double big = 0.0;
        int irow = -1, icol = -1;

        // 1. 寻找主元 (Pivot)
        for (j = 0; j < n; j++) {
            if (ipiv[j] != 1) {
                for (k = 0; k < n; k++) {
                    if (ipiv[k] == 0) {
                        // 行主序索引: A[j][k] -> A[j*n + k]
                        if (fabs(A[j * n + k]) >= big) {
                            big = fabs(A[j * n + k]);
                            irow = j;
                            icol = k;
                        }
                    }
                }
            }
        }

        if (irow == -1 || icol == -1 || big <= scale * 1e-12) return -1; // 奇异矩阵
        ipiv[icol]++;

        // 2. 交换行 (将主元行 irow 换到 icol 位置)
        // 这样主元就位于 A[icol][icol]
        if (irow != icol) {
            for (int l = 0; l < n; l++) {
                swap(A[irow * n + l], A[icol * n + l]);
            }
        }

        indxr[i] = irow;
        indxc[i] = icol;

        // 检查对角线元素
        if (A[icol * n + icol] == 0.0) return -1;

        // 3. 归一化当前行 (Pivot Row)
        double pivinv = 1.0 / A[icol * n + icol];
        A[icol * n + icol] = 1.0;
        for (int l = 0; l < n; l++) A[icol * n + l] *= pivinv;

        // 4. 消元 (Elimination)
        for (int ll = 0; ll < n; ll++) {
            if (ll != icol) { // 跳过当前主元行
                double dum = A[ll * n + icol];
                A[ll * n + icol] = 0.0;
                for (int l = 0; l < n; l++) A[ll * n + l] -= A[icol * n + l] * dum;
            }
        }
    }

    // 5. 恢复列交换 (因为Gauss-Jordan全主元法交换了列，最后需要换回来)
    for (int l = n - 1; l >= 0; l--) {
        if (indxr[l] != indxc[l]) {
            for (int k = 0; k < n; k++) {
                // 交换第 k 行的两个元素: A[k][indxr[l]] 和 A[k][indxc[l]]
                swap(A[k * n + indxr[l]], A[k * n + indxc[l]]);
            }
        }
    }

    return 0;
}

bool leastSquares(const double* H, const double* residual, const double* weights,
                  int rows, int cols, double* solution) {
    if (cols<1 || rows<cols) return false;
    std::vector<double> a(rows*cols), b(rows);
    double scale=0;
    for (int i=0; i<rows; ++i) {
        if (!std::isfinite(weights[i]) || weights[i]<=0 || !std::isfinite(residual[i])) return false;
        const double w=std::sqrt(weights[i]); b[i]=residual[i]*w;
        for (int j=0; j<cols; ++j) {
            a[i*cols+j]=H[i*cols+j]*w;
            if (!std::isfinite(a[i*cols+j])) return false;
            scale=std::max(scale,std::fabs(a[i*cols+j]));
        }
    }
    for (int j=0; j<cols; ++j) {
        double norm=0;
        for (int i=j; i<rows; ++i) norm=std::hypot(norm,a[i*cols+j]);
        if (norm<=scale*1e-10 || !std::isfinite(norm)) return false;
        std::vector<double> u(rows-j);
        for (int i=j; i<rows; ++i) u[i-j]=a[i*cols+j];
        u[0]+=std::copysign(norm,u[0]);
        double unorm=0;
        for (double x:u) unorm=std::hypot(unorm,x);
        for (double& x:u) x/=unorm;
        for (int k=j; k<cols; ++k) {
            double dot=0;
            for (int i=j; i<rows; ++i) dot+=u[i-j]*a[i*cols+k];
            for (int i=j; i<rows; ++i) a[i*cols+k]-=2*u[i-j]*dot;
        }
        double dot=0;
        for (int i=j; i<rows; ++i) dot+=u[i-j]*b[i];
        for (int i=j; i<rows; ++i) b[i]-=2*u[i-j]*dot;
    }
    for (int j=cols-1; j>=0; --j) {
        double v=b[j];
        for (int k=j+1; k<cols; ++k) v-=a[j*cols+k]*solution[k];
        solution[j]=v/a[j*cols+j];
        if (!std::isfinite(solution[j])) return false;
    }
    return true;
}
