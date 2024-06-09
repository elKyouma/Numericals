#include <gtest/gtest.h>
#include <valarray>
#include "MatrixSolver.h"
#include "PivotingStrategy.h"
#include "vector.h"
#include "MatrixDecomposer.h"

using namespace numericals;

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

TEST(MatrixEquationSolver, SolveHighTriangular)
{
    matrix<real> A{3, 3, {1.0, 2.0, 3.0, 0.0, 5.0, 6.0, 0.0, 0.0, 9.0}};
    vector<real> b{14.0, 28.0, 27.0};

    auto x = solve_high_trian_matrix_eq(A, b);

    ASSERT_EQ(x.GetSize(), 3);
    expect_valarray_equals<real>(x, std::valarray<real>{1.0, 2.0, 3.0});
}

TEST(MatrixEquationSolver, SolveLowTriangular)
{
    matrix<real> A{3, 3, {  1.0, 0.0, 0.0, 
                            2.0, 3.0, 0.0, 
                            4.0, 5.0, 6.0}};
    vector<real> b{1.0, 2.0, 10.0};

    auto x = solve_low_trian_matrix_eq(A, b);

    ASSERT_EQ(x.GetSize(), 3);
    expect_valarray_equals<real>(x, std::valarray<real>{1.0, 0.0, 1.0});
}

TEST(MatrixEquationSolver, SolveGauss)
{
    matrix<real> A{3, 3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 9.0, 8.0, 8.0}};
    vector<real> b{14.0, 32.0, 50.0};

    auto x = solve_matrix_eq_gauss(A, b);

    ASSERT_EQ(x.GetSize(), 3);
    expect_valarray_equals<real>(x, std::valarray<real>{2.0, 0.0, 4.0});
}

TEST(MatrixEquationSolver, SolveJordan)
{
    matrix<real> A{3, 3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 9.0, 8.0, 8.0}};
    vector<real> b{14.0, 32.0, 50.0};

    auto x = solve_matrix_eq_jordan(A, b);

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
    auto x1 = solve_matrix_eq_gauss(A, b, PartialPivotingStragegy());
    auto x2 = solve_matrix_eq_jordan(A, b, PartialPivotingStragegy());
    
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
    auto x1 = solve_matrix_eq_gauss(A, b, FullPivotingStragegy());
    auto x2 = solve_matrix_eq_jordan(A, b, FullPivotingStragegy());
    
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

TEST(MatrixEquationSolver, LDLT_Decomposition)
{
    matrix<real> A{3, 3, {  1.0, 0.0, 3.0,
                            0.0, 6.0, 6.0,
                            3.0, 6.0, 5.0 }};
    matrix<real> expected_lu{3, 3, {1.0, 0.0, 3.0,
                                    0.0, 6.0, 1.0,
                                    3.0, 1.0, -10.0 }};
    A = ldlt_decomposition(A);

    expect_matrix_equals(A, expected_lu); 
}

TEST(MatrixEquationSolver, LLT_Decomposition)
{
    matrix<real> A{3, 3, {  4.0, 12.0, -16.0,
                            12.0, 37.0, -43.0,
                            -16.0, -43.0, 98.0 }};
    
    matrix<real> expected_lu{3, 3, {2.0, 6.0, -8.0,
                                    6.0, 1.0, 5.0,
                                    -8.0, 5.0, 3.0 }};
    A = llt_decomposition(A);

    expect_matrix_equals(A, expected_lu); 
}

TEST(MatrixEquationSolver, solveOverdeterminedMatrix)
{
    matrix<real> A{3, 4, {  1.0, 0.0, 0.0,
                            0.0, 2.0, 0.0,
                            2.0, 0.0, 1.0,
                            0.0, 0.0, 1.0}};
    vector<real> b{1.0, 2.0, 3.0, 0.0}; 
    //vector<real> b{72.0, 0.0, 288.0};
    //vector<real> expected{4162.0, -1136.0, 184.0};
    auto x = solve_overdetermined_matrix(A, b, solve_matrix_eq_jordan);
    std::cout <<x;
    //expect_valarray_equals((std::valarray<real>)x, (std::valarray<real>)expected); 
}

TEST(MatrixEquationSolver, solveLLT_Matrix)
{
    matrix<real> A{3, 3, {  4.0, 12.0, -16.0,
                            12.0, 37.0, -43.0,
                            -16.0, -43.0, 98.0 }};
    
    vector<real> b{72.0, 0.0, 288.0};
    vector<real> expected{4162.0, -1136.0, 184.0};
    auto x = solve_matrix_eq_with_llt_decomposition(A,b);

    expect_valarray_equals((std::valarray<real>)x, (std::valarray<real>)expected); 
}

TEST(MatrixEquationSolver, SolveLU_Matrix)
{
    matrix<real> A{3, 3, {  1.0, 0.0, 3.0,
                            0.0, 6.0, 6.0,
                            3.0, 6.0, 5.0 }};
    vector<real> b{1.0, 0.0, 13.0};
    vector<real> expected{4.0, 1.0, -1.0};
    vector<real> x = solve_matrix_eq_with_lu_decomposition(A, b);

    expect_valarray_equals((std::valarray<real>)x, (std::valarray<real>)expected); 
}

TEST(MatrixEquationSolver, SolveLDLT_Matrix)
{
    matrix<real> A{3, 3, {  1.0, 0.0, 3.0,
                            0.0, 6.0, 6.0,
                            3.0, 6.0, 5.0 }};
    
    vector<real> b{1.0, 0.0, 13.0};
    vector<real> x = solve_matrix_eq_with_ldlt_decomposition(A, b);
    vector<real> expected{4.0, 1.0, -1.0};
    expect_valarray_equals((std::valarray<real>)x, (std::valarray<real>)expected); 
}

TEST(MatrixEquationSolver, SolveTridiagonalMatrix)
{
    vector<real> a1 { 1.0, 6.0, 2.0 };
    vector<real> a2 { 1.0, 4.0, 6.0, 2.0 };
    vector<real> a3 { 2.0, 4.0, 8.0 };
    vector<real> b {14.0, 0.0, 14.0, 14.0};

    vector<real> expected{28.0, -7.0, 0.0, 7.0};
    auto result = solve_tridiagonal_matrix_eq({a1,a2,a3}, b);
    expect_valarray_equals((std::valarray<real>)result, (std::valarray<real>)expected); 
}

