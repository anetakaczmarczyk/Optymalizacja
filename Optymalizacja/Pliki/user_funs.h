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
matrix df2(double, matrix, matrix = NAN, matrix = NAN);
matrix ff2T(matrix, matrix = NAN, matrix = NAN);
matrix ff2R(matrix, matrix = NAN, matrix = NAN);

//lab 3
matrix df3(double, matrix, matrix = NAN, matrix = NAN);
matrix ff3T_out(matrix, matrix = NAN, matrix = NAN);
matrix ff3T_in(matrix, matrix = NAN, matrix = NAN);
matrix ff3R(matrix, matrix = NAN, matrix = NAN);

//lab 4
matrix ff4T(matrix x, matrix ud1= NAN, matrix ud2= NAN);
matrix gfT(matrix x, matrix ud1= NAN, matrix ud2= NAN);
matrix hfT(matrix x, matrix ud1= NAN, matrix ud2= NAN);
matrix ff4R(matrix theta, matrix X= NAN, matrix Y= NAN);
matrix gf4R(matrix theta, matrix X= NAN, matrix Y= NAN);
double sygmoid(matrix theta, matrix x);