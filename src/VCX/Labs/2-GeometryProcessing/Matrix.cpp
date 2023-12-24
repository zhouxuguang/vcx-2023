//
// Created by zhouxuguang on 2023/12/2.
//

#include "Matrix.h"
#include <math.h>
#include <algorithm>
#include <limits>

Matrix::Matrix(int width, int height)
{
    if (!width || !height)
    {
        return;
    }

    mData = new double[width * height];
    memset(mData, 0, width * height * sizeof(double));
    mWidth = width;
    mHeight = height;
}

Matrix::~Matrix()
{
    if (mData)
    {
        delete []mData;
        mData = nullptr;
    }
}

double *Matrix::operator[](int row)
{
    return mData + (row * mWidth);
}

double *Matrix::operator[](int row) const
{
    return mData + (row * mWidth);
}

void Matrix::SwapRows(int row1, int row2)
{
    double *row1Data = mData + (row1 * mWidth);
    double *row2Data = mData + (row2 * mWidth);
    for (int i = 0; i < mWidth; i++)
    {
        std::swap(row1Data[i], row2Data[i]);
    }
}

int Matrix::GetColumnMaxAbsElementIndex(int col, int startRow, int endRow)
{
    int rowIndex = startRow;
    double maxValue = std::numeric_limits<double>::min();
    for (int i = startRow; i <= endRow; ++i)
    {
        double value = *(mData + (i * mWidth + col));
        if (fabs(value) >= maxValue)
        {
            maxValue = value;
            rowIndex = i;
        }
    }
    return rowIndex;
}

void Matrix::AddKScaleToOther(int row1, double scaleFactor, int row2)
{

}

