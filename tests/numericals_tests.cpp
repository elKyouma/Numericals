#include <gtest/gtest.h>
#include "numericals.h"
#include "vector.h"

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

TEST(MatrixEquationSolver, Solve)
{
    matrix<real> A{3, 3, {1.0, 2.0, 3.0, 0.0, 5.0, 6.0, 0.0, 0.0, 9.0}};
    vector<real> b{14.0, 28.0, 27.0};

    auto x = solve_triangular_matrix_equation(A, b);

    ASSERT_EQ(x.size(), 3);
    ASSERT_FLOAT_EQ(x[0], 1.0);
    ASSERT_FLOAT_EQ(x[1], 2.0);
    ASSERT_FLOAT_EQ(x[2], 3.0);
}
