#pragma once
#include "matrix.h"
#include "numerical_types.h"
#include "PivotingStrategy.h"

namespace numericals {

matrix<real> lu_decomposition(matrix<real> a, PivotingStrategy&& strategy = NoPivotingStragegy());
matrix<real> ldlt_decomposition(const matrix<real>& a, PivotingStrategy&& strategy = NoPivotingStragegy());
matrix<real> ldl_decomposition(const matrix<real>& a, PivotingStrategy&& strategy = NoPivotingStragegy());

}
