#include <gtest/gtest.h>
#include "utils.h"
#include "matrix.h"

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
