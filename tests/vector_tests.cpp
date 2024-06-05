#include "vector.h"
#include <gtest/gtest.h>

TEST(Vectors, create)
{
    std::array<double, 4> init {2.0, 3.0, 8.0, 5.0};
    vector<double> vec{2.0, 3.0, 8.0, 5.0};

    for(size_t i = 0; i < vec.size(); i++)
        EXPECT_EQ(vec[i], init[i]);
}
