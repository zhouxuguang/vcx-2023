//
// Created by zhouxuguang on 2023/12/2.
//

#include "LinearSystem.h"
#include <limits>
#include <iostream>
#include <assert.h>

void LinearSystem::GaussElimination(const Matrix &A, const Vector &B, Vector &X)
{
    //1 构造增广矩阵
    Matrix extendMat = Matrix(A.GetWith() + 1, A.GetHeight());
    for (int i = 0; i < A.GetHeight(); ++i)
    {
        memcpy(extendMat[i], A[i], A.GetWith() * sizeof(double));
        extendMat[i][A.GetWith()] = B[i];
    }

    int height = A.GetHeight();
    int width = extendMat.GetWith();

    //每一行逐步进行消元，即遍历对角元
    for (int k = 0; k < height - 1; ++k)
    {
        int maxIndex = extendMat.GetColumnMaxAbsElementIndex(k, k, height - 1);

        //当前行的第一个元素，akk不是绝对值最大元素，就需要将当前元素所在行和绝对值最大的行进行交换
        if (maxIndex != k)
        {
            extendMat.SwapRows(k, maxIndex);
        }

        //对第k行下方的所有行进行遍历
        for (int i = k + 1; i < height; ++i)
        {
            //计算好每一行的需要乘以的系数
            extendMat[i][k] = extendMat[i][k] / extendMat[k][k];

            //对第i行，从k+1列开始遍历
            for (int j = k + 1; j < width; ++j)
            {
                extendMat[i][j] = extendMat[i][j] - extendMat[i][k] * extendMat[k][j];
            }
        }
    }

    //开始回代求解
    for (int i = height - 1; i >= 0; --i)
    {
        double sum = 0;
        for (int j = i + 1; j < height; ++j)
        {
            sum += extendMat[i][j] * X[j];
        }

        X[i] = (extendMat[i][width - 1] - sum) / extendMat[i][i];
    }
}

void LinearSystem::JacobiIteration(const Matrix &A, const Vector &B, Vector &X)
{
    int width = A.GetWith();
    int height = A.GetHeight();
    assert(width == height);

    Vector X0;
    X0.resize(X.size());
    memset(X0.data(), 0, X0.size() * sizeof(double));

    Vector diffVec;
    diffVec.resize(X.size());

    while (true)
    {
        for (int i = 0; i < width; ++i)
        {
            double sum = 0.0;
            for (int j = 0; j < width; ++j)
            {
                if (i == j)
                {
                    continue;
                }
                sum += A[i][j] * X0[j];
            }

            // xi(k1)a ( ax(k) ax(k)b)
            X[i] = (B[i] - sum) / A[i][i];
        }

        //判断X0和X相差的阈值小于一定量，停止迭代
        double sum = 0.0;
        for (int i = 0; i < width; ++i)
        {
            diffVec[i] = X[i] - X0[i];
            sum += diffVec[i] * diffVec[i];
        }

        std::cout << "x1 = " << X[0] << " x2 = " << X[1] << " x3 = " << X[2] << std::endl;

        if (sum <= std::numeric_limits<double>::epsilon())
        {
            break;
        }

        //更新X0
        memcpy(X0.data(), X.data(), X.size() * sizeof(double));
    }
}

void LinearSystem::GaussSeidelIteration(const Matrix &A, const Vector &B, Vector &X, int iterCount)
{
    int width = A.GetWith();
    int height = A.GetHeight();
    assert(width == height);

    Vector X0;
    X0.resize(X.size());
    memset(X0.data(), 0, X0.size() * sizeof(double));

    Vector diffVec;
    diffVec.resize(X.size());
    
    int curIterCount = 0;

    while (true)
    {
        curIterCount ++;
        for (int i = 0; i < width; ++i)
        {
            double sum = 0.0;

            //在计算xi的时候，新的xi前面的都计算好了，直接用新的来计算
            for (int j = 0; j < i; ++j)
            {
                sum += A[i][j] * X[j];
            }

            for (int j = i + 1; j < width; ++j)
            {
                sum += A[i][j] * X0[j];
            }

            X[i] = (B[i] - sum) / A[i][i];
        }

        //判断X0和X相差的阈值小于一定量，停止迭代
        double sum = 0.0;
        for (int i = 0; i < width; ++i)
        {
            diffVec[i] = X[i] - X0[i];
            sum += diffVec[i] * diffVec[i];
        }

        std::cout << "x1 = " << X[0] << " x2 = " << X[1] << " x3 = " << X[2] << std::endl;

        if (sum <= std::numeric_limits<double>::epsilon() || curIterCount >= iterCount)
        {
            break;
        }

        //更新X0
        memcpy(X0.data(), X.data(), X.size() * sizeof(double));
    }
    
    printf("gauss -seidel iter count = %d\n", iterCount);
}
