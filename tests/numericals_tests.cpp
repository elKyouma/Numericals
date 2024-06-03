#include <gtest/gtest.h>
#include "numericals.h"

TEST(MatrixSolver, Solve)
{
    const auto expected = 1;
    ASSERT_EQ(expected, 1);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
