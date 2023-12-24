//
// Created by zhouxuguang on 2023/12/2.
//

#ifndef NUMERICALANALYSIS_MATRIX_H
#define NUMERICALANALYSIS_MATRIX_H

//矩阵类
class Matrix
{
public:
    Matrix(int width, int height);

    ~Matrix();

    int GetWith() const
    {
        return mWidth;
    }

    int GetHeight() const
    {
        return mHeight;
    }

    double * operator [] (int row);

    double * operator [] (int row) const;

    /**
     * 交换两行的数据
     * @param row1
     * @param row2
     */
    void SwapRows(int row1, int row2);

    /**
     * 将第row1行的scaleFactor加到第row2行，这个先不用吧
     * @param row1
     * @param scaleFactor
     * @param row2
     */
    void AddKScaleToOther(int row1, double scaleFactor, int row2);

    /**
     * 获得第col的元素最大值所在的位置
     * @param col
     * @param startRow
     * @param endRow
     * @return
     */
    int GetColumnMaxAbsElementIndex(int col, int startRow, int endRow);
private:
    int mWidth = 0;   //宽，列数
    int mHeight = 0;  //高，行数
    double *mData = nullptr;

};


#endif //NUMERICALANALYSIS_MATRIX_H
