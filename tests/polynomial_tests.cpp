#include "PolynomialSolver.h"
#include <cmath>
#include <gtest/gtest.h>
#include <numeric>

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

TEST(Polynomial, GetCzebyszewZeros)
{
    vector<real> expected{-4.82963, -3.53553, -1.2941, 1.2941, 3.53553, 4.82963};
    auto vec1 = get_chebyshev_polynomial_zeros(6, -5, 5);

    double abs_error = 0.00001;
    EXPECT_NEAR(expected[0], vec1[0], abs_error);
    EXPECT_NEAR(expected[1], vec1[1], abs_error);
    EXPECT_NEAR(expected[2], vec1[2], abs_error);
    EXPECT_NEAR(expected[3], vec1[3], abs_error);
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

TEST(Approximation, LeastSquareApprox_Simple)
{
    std::vector<real> x (5);
    std::iota(x.begin(), x.end(), -2);
    std::vector<real> y {4, 1, 0, 1, 4};
    vector<real> expected{0, 0, 1, 0};

    double abs_error = 0.0001;
    
    auto polynomial = get_polynomial_approximation(x, y, 3);
    EXPECT_NEAR(polynomial[0], expected[0], abs_error);
    EXPECT_NEAR(polynomial[1], expected[1], abs_error);
    EXPECT_NEAR(polynomial[2], expected[2], abs_error);
    EXPECT_NEAR(polynomial[3], expected[3], abs_error);

}

TEST(Approximation, LeastSquareApprox_Advanced)
{
    std::vector<real> x (5);
    std::iota(x.begin(), x.end(), -3);
    std::vector<real> y {4, 1, 0, 1, 4};
    vector<real> expected{1, 2, 1, 0};

    double abs_error = 0.0001;
    
    auto polynomial = get_polynomial_approximation(x, y, 3);
    auto polynomial2 = get_polynomial_approximation(x, y, 3);
    EXPECT_NEAR(polynomial[0], expected[0], abs_error);
    EXPECT_NEAR(polynomial[1], expected[1], abs_error);
    EXPECT_NEAR(polynomial[2], expected[2], abs_error);
    EXPECT_NEAR(polynomial[3], expected[3], abs_error);

}

TEST(Interpolation, LagrangeInterpolation)
{
    std::vector<real> x (5);
    std::iota(x.begin(), x.end(), -3);
    std::vector<real> y {4, 1, 0, 1, 4};
    double abs_error = 0.0001;
    
    auto func = get_lagrange_interpolation(x, y);

    EXPECT_NEAR(func(-3), 4.0, abs_error);
    EXPECT_NEAR(func(-2), 1.0, abs_error);
    EXPECT_NEAR(func(-1), 0.0, abs_error);
    EXPECT_NEAR(func(0), 1.0, abs_error);
    EXPECT_NEAR(func(1), 4.0, abs_error);
}

TEST(Interpolation, NewtonInterpolation)
{
    std::vector<real> x (5);
    std::iota(x.begin(), x.end(), -3);
    std::vector<real> y {4, 1, 0, 1, 4};
    double abs_error = 0.0001;
    
    auto func = get_newton_interpolation(x, y);

    EXPECT_NEAR(func(-3), 4.0, abs_error);
    EXPECT_NEAR(func(-2), 1.0, abs_error);
    EXPECT_NEAR(func(-1), 0.0, abs_error);
    EXPECT_NEAR(func(0), 1.0, abs_error);
    EXPECT_NEAR(func(1), 4.0, abs_error);
}

