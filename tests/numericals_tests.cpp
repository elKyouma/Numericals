#include <gtest/gtest.h>
#include <valarray>
#include "numericals.h"
#include "vector.h"


template <typename T> requires std::is_arithmetic_v<T>
void expect_valarray_equals(std::valarray<T> ar1, std::valarray<T> ar2)
{
    for(size_t i = 0; i < ar1.size(); i++)
        EXPECT_NEAR(ar1[i], ar2[i], 10e-7);
}

template <typename T> requires std::is_arithmetic_v<T>
void expect_matrix_equals(matrix<T> a1, matrix<T> a2)
{
    for(size_t i = 0; i < a1.GetSizeY() * a1.GetSizeX(); i++)
        EXPECT_NEAR(a1.GetElement(i), a2.GetElement(i), 10e-7);
}

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

TEST(MatrixEquationSolver, SolveTriangular)
{
    matrix<real> A{3, 3, {1.0, 2.0, 3.0, 0.0, 5.0, 6.0, 0.0, 0.0, 9.0}};
    vector<real> b{14.0, 28.0, 27.0};

    auto x = solve_triangular_matrix_equation(A, b);

    ASSERT_EQ(x.GetSize(), 3);
    expect_valarray_equals<real>(x, std::valarray<real>{1.0, 2.0, 3.0});
}

TEST(MatrixEquationSolver, SolveGauss)
{
    matrix<real> A{3, 3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 9.0, 8.0, 8.0}};
    vector<real> b{14.0, 32.0, 50.0};

    auto x = solve_matrix_equation_gauss(A, b, MatrixFlag::NORMAL);

    ASSERT_EQ(x.GetSize(), 3);
    expect_valarray_equals<real>(x, std::valarray<real>{2.0, 0.0, 4.0});
}

TEST(MatrixEquationSolver, SolveJordan)
{
    matrix<real> A{3, 3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 9.0, 8.0, 8.0}};
    vector<real> b{14.0, 32.0, 50.0};

    auto x = solve_matrix_equation_jordan(A, b, MatrixFlag::NORMAL);

    ASSERT_EQ(x.GetSize(), 3);
    expect_valarray_equals<real>(x, std::valarray<real>{2.0, 0.0, 4.0});
}

TEST(MatrixEquationSolver, PartialSelection)
{
    matrix<real> A{3, 3, {  1.0, 9.0, 2.0,
                            2.0, 5.0, 7.0,
                            3.0, 8.0, 3.0}};
    
    vector<real> b {1.0, 2.0, 3.0};
    vector<real> solution {1.0, 0.0, 0.0};
    auto x1 = solve_matrix_equation_gauss(A, b, MatrixFlag::PARTIAL_SELECT);
    auto x2 = solve_matrix_equation_jordan(A, b, MatrixFlag::PARTIAL_SELECT);
    
    expect_valarray_equals<real>(x1, solution);
    expect_valarray_equals<real>(x2, solution);
}

TEST(MatrixEquationSolver, FullSelection)
{
    matrix<real> A{3, 3, {  8.0, 8.0, 8.0,
                            4.0, 5.0, 7.0,
                            1.0, 2.0, 3.0 }};
    
    vector<real> b {8.0, 1.0, 0.0};
    vector<real> solution {0.0, 3.0, -2.0};
    auto x1 = solve_matrix_equation_gauss(A, b, MatrixFlag::FULL_SELECT);
    auto x2 = solve_matrix_equation_jordan(A, b, MatrixFlag::FULL_SELECT);
    
    expect_valarray_equals<real>(x1, solution);
    expect_valarray_equals<real>(x2, solution);
}

TEST(MatrixEquationSolver, LU_Decomposition)
{
    matrix<real> A{3, 3, {  1.0, 2.0, 3.0,
                            4.0, 5.0, 6.0,
                            7.0, 8.0, 9.0 }};
    matrix<real> expected_lu{3, 3, {1.0, 2.0, 3.0,
                                    4.0, -3.0, -6.0,
                                    7.0, 2.0, 0.0 }};
    A = lu_decomposition(A);

    expect_matrix_equals(A, expected_lu); 
}

