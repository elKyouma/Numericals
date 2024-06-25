#pragma once
#include "matrix.h"
#include "numerical_types.h"
#include "PivotingStrategy.h"
#include <utility>

namespace numericals {

matrix<real> lu_decomposition(matrix<real> a, PivotingStrategy&& strategy = NoPivotingStragegy());
matrix<real> ldlt_decomposition(const matrix<real>& a);
matrix<real> llt_decomposition(const matrix<real>& a);
std::pair<matrix<real>, matrix<real>> qr_decomposition(const matrix<real>& a);

}
