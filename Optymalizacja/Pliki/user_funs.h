#pragma once

#include"ode_solver.h"


matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);

// lab 1
matrix ff1(matrix, matrix = NAN, matrix = NAN);
matrix f1R(matrix, matrix = NAN, matrix = NAN);
matrix df1(double, matrix, matrix = NAN, matrix = NAN);

// lab 2
matrix ff2T(matrix, matrix, matrix = NAN);