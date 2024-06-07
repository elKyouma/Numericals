#include <gtest/gtest.h>
#include <valarray>
#include "utils.h"
#include "matrix.h"

TEST(Utils, SwapValArrays)
{
    std::valarray<double> x{1.0, 2.0, 3.0};
    std::valarray<double> y{4.0, 5.0, 6.0};

    swap_slices(x, y);
    EXPECT_EQ(x[0], 4);
    EXPECT_EQ(x[1], 5);
    EXPECT_EQ(x[2], 6);

    EXPECT_EQ(y[0], 1);
    EXPECT_EQ(y[1], 2);
    EXPECT_EQ(y[2], 3);
}

TEST(Utils, SwapSliceArrays)
{
    std::valarray<double> x{1.0, 2.0, 3.0};
    std::valarray<double> y{4.0, 5.0, 6.0};
    std::slice_array<double>xx = x[std::slice(0, 3, 1)];
    std::slice_array<double>yy = y[std::slice(0, 3, 1)];

    swap_slices(xx, yy);
    EXPECT_EQ(x[0], 4);
    EXPECT_EQ(x[1], 5);
    EXPECT_EQ(x[2], 6);

    EXPECT_EQ(y[0], 1);
    EXPECT_EQ(y[1], 2);
    EXPECT_EQ(y[2], 3);
}

TEST(Utils, FindIndexOfMaxElementInTheValarray)
{
    matrix<double> A{3, 3, {1.0, 2.0, 3.0,
                            2.0, 1.0, 4.0,
                            3.0, 0.0, 1.0}};
    
    EXPECT_EQ(2, find_index_of_column_max(A.GetColumn(0), 0, 3));
    EXPECT_EQ(2, find_index_of_column_max(A.GetColumn(0), 1, 3));
    EXPECT_EQ(2, find_index_of_column_max(A.GetColumn(0), 2, 3));

    EXPECT_EQ(0, find_index_of_column_max(A.GetColumn(1), 0, 3));
    EXPECT_EQ(1, find_index_of_column_max(A.GetColumn(1), 1, 3));
    EXPECT_EQ(2, find_index_of_column_max(A.GetColumn(1), 2, 3));

    EXPECT_EQ(1, find_index_of_column_max(A.GetColumn(2), 0, 3));
    EXPECT_EQ(1, find_index_of_column_max(A.GetColumn(2), 1, 3));
    EXPECT_EQ(2, find_index_of_column_max(A.GetColumn(2), 2, 3));
}
