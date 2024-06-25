#include "PolynomialSolver.h"
#include <cmath>
#include <gtest/gtest.h>

TEST(Polynomials, Solve)
{
    //x^2 + 2x + 1
    real a[3];
    a[0] = 1;
    a[1] = 2;
    a[2] = 1;
    
    auto result = solve_polynomial(a, 3.0); 
    ASSERT_FLOAT_EQ(result, 16.0);
}

TEST(Polynomials, SolveHorner)
{
    //x^2 + 2x + 1
    real a[3];
    a[0] = 1;
    a[1] = 2;
    a[2] = 1;
    
    auto result = solve_polynomial_horner(a, 3.0); 
    ASSERT_FLOAT_EQ(result, 16.0);
}

TEST(Functions, FindZero)
{
    auto func = [](double x){return exp(x) - 1;};
    EXPECT_FLOAT_EQ(0.f, find_function_zero_with_bisection( func, -1, 1));
}
