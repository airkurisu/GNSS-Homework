/*
* include/matrix.h - 矩阵运算库声明
* 项目: MYGNSS10
* 说明: 采用行主序 (Row-Major) 存储，兼容 C++ 标准数组布局
*/

#ifndef MYGNSS10_MATRIX_H
#define MYGNSS10_MATRIX_H

/**
 * 矩阵乘法: C = alpha * op(A) * op(B) + beta * C
 * @param tr     转置标志 ("NN": A不转B不转, "NT": A不转B转, "TN": A转B不转, "TT": A转B转)
 * @param n      op(A)的行数，即 C 的行数
 * @param k      op(A)的列数，即 op(B) 的行数 (公共维度)
 * @param m      op(B)的列数，即 C 的列数
 * @param alpha  乘法系数
 * @param A      矩阵 A (一维数组表示)
 * @param B      矩阵 B (一维数组表示)
 * @param beta   加法系数 (通常为 0.0)
 * @param C      结果矩阵 C (n x m)
 */
void matmul(const char *tr, int n, int k, int m, double alpha,
            const double *A, const double *B, double beta, double *C);

/**
 * 矩阵求逆 (高斯-约旦消元法)
 * 注意：直接在原矩阵上操作，计算后 A 变为其逆矩阵
 * @param A      输入矩阵 (n x n)，输出时变为逆矩阵
 * @param n      矩阵维数
 * @return       0: 成功, -1: 奇异矩阵(不可逆)
 */
int matinv(double *A, int n);

#endif // MYGNSS10_MATRIX_H