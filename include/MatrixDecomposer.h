#pragma once
#include "matrix.h"
#include "numerical_types.h"

matrix<real> lu_decomposition(matrix<real> a, MatrixFlag flag = NORMAL);
matrix<real> ldlt_decomposition(const matrix<real>& a, MatrixFlag flag = NORMAL);
matrix<real> ldl_decomposition(const matrix<real>& a, MatrixFlag flag = NORMAL);
