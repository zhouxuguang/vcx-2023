//
// Created by zhouxuguang on 2023/12/2.
//

#ifndef NUMERICALANALYSIS_LINEARSYSTEM_H
#define NUMERICALANALYSIS_LINEARSYSTEM_H

#include "PreDefine.h"
#include "Matrix.h"

class LinearSystem
{
public:
    /**
     * 列主元高斯消元法
     * @param A
     * @param B
     * @param X
     */
    static void GaussElimination(const Matrix& A, const Vector& B, Vector &X);

    /**
     * 雅可比迭代解线性方程组
     * @param A
     * @param B
     * @param X
     */
    static void JacobiIteration(const Matrix& A, const Vector& B, Vector &X);

    /**
     * 高斯塞德尔迭代求解线性方程组
     * @param A
     * @param B
     * @param X
     */
    static void GaussSeidelIteration(const Matrix& A, const Vector& B, Vector &X);
};


#endif //NUMERICALANALYSIS_LINEARSYSTEM_H
