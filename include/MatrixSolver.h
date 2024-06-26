#pragma once
#include "numerical_types.h"
#include "matrix.h"
#include "vector.h"
#include "PivotingStrategy.h"
#include <functional>

namespace numericals {
//w = x^n-1*a_n + x^n-2*a_n-1 + a_0 
vector<real> solve_high_trian_matrix_eq(const matrix<real>& a, const vector<real>& b, bool assumeDiagonalOnes = false);
vector<real> solve_low_trian_matrix_eq(const matrix<real>& a, const vector<real>& b, bool assumeDiagonalOnes = false);
vector<real> solve_overdetermined_matrix(const matrix<real>& a, const vector<real>& b, std::function<vector<real>(matrix<real>, vector<real>, PivotingStrategy&&)> algorithm, PivotingStrategy&& strategy = NoPivotingStragegy());

vector<real> solve_matrix_eq_gauss( matrix<real> a, vector<real> b, PivotingStrategy&& strategy = NoPivotingStragegy());
vector<real> solve_matrix_eq_jordan( matrix<real> a, vector<real> b,  PivotingStrategy&& strategy = NoPivotingStragegy());
vector<real> solve_tridiagonal_matrix_eq( std::array<vector<real>, 3> a, vector<real> b);
vector<real> solve_matrix_eq_with_qr_decomposition(const matrix<real>& a, const vector<real>& b);
vector<real> solve_matrix_eq_with_lu_decomposition(const matrix<real>& a, const vector<real>& b, PivotingStrategy&& strategy = NoPivotingStragegy());
vector<real> solve_matrix_eq_with_ldlt_decomposition(const matrix<real>& a, const vector<real>& b);
vector<real> solve_matrix_eq_with_llt_decomposition(const matrix<real>& a, const vector<real>& b);

}
