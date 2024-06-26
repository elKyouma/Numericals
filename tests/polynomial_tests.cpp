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
    double abs_error = 0.0001;

    auto func1 = [](double x){return exp(x) - 1;};
    EXPECT_NEAR(0.f, find_function_zero_with_bisection( func1, -8, 1.2), abs_error);
    EXPECT_NEAR(0.f, find_function_zero_with_falsi( func1, -8, 1.2), abs_error);
    
    auto func2 = [](double x){return (x - 5.)*(x + 30.);};
    EXPECT_NEAR(5.f, find_function_zero_with_bisection( func2, -1, 9.2), abs_error);
    EXPECT_NEAR(5.f, find_function_zero_with_falsi( func2, -1, 9.2), abs_error);
    EXPECT_NEAR(5.f, find_function_zero_with_secant( func2, -1, 9.2), abs_error);
    EXPECT_NEAR(5.f, find_function_zero_with_newton_raphson( func2, -1, 9.2 ), abs_error);
}

