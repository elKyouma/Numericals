#pragma once

#include <stack>
#include <utility>

typedef float real;
//typedef double real;

typedef std::stack<std::pair<size_t, size_t>> permutation_stack;

enum MatrixFlag
{
    NORMAL,
    PARTIAL_SELECT,
    FULL_SELECT
};

