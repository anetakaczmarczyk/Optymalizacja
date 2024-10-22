#pragma once

#include"ode_solver.h"

matrix ff1(matrix, matrix = NAN, matrix = NAN);
matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);

