#include "matrix.h"
#include <gtest/gtest.h>

TEST(Matrices, create)
{
    std::array<double, 2> col1 {2.0, 8.0};
    std::array<double, 2> col2 {3.0, 5.0};
    std::array<double, 2> row1 {2.0, 3.0};
    std::array<double, 2> row2 {8.0, 5.0};

    matrix<double> mat{2, 2, {2.0, 3.0, 8.0, 5.0}};
    EXPECT_EQ(mat.GetRow(0)[0], row1[0]);
    EXPECT_EQ(mat.GetRow(0)[1], row1[1]);
    EXPECT_EQ(mat.GetRow(1)[0], row2[0]);
    EXPECT_EQ(mat.GetRow(1)[1], row2[1]);

    EXPECT_EQ(mat.GetColumn(0)[0], col1[0]);
    EXPECT_EQ(mat.GetColumn(0)[1], col1[1]);
    EXPECT_EQ(mat.GetColumn(1)[0], col2[0]);
    EXPECT_EQ(mat.GetColumn(1)[1], col2[1]);
}
