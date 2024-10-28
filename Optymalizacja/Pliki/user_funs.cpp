#include"user_funs.h"

matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
	int n = get_len(Y[0]);
	double teta_max = Y[1](0, 0);
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);
	double m = 1, l = 0.5, b = 0.5, g = 9.81;
	double I = m*pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1))*ud2(0) - m*g*l*sin(Y(0)) - b*Y(1)) / I;
	return dY;
}

// tu implementowac funkcje testową
matrix ff1(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = -cos(0.1 * m2d(x)) * exp (-1.0 * pow((0.1 * m2d(x) - 2 * M_PI), 2)) + 0.002 * pow(0.1 * m2d(x), 2);
	return y;
}
matrix f1R(matrix x , matrix ud1, matrix ud2){
    // //dla tego nie uzywamy metody ekspansji
    // matrix y;
    // matrix Y0=matrix(3,new double[]{5,1,20});
    // matrix* Y= solve_ode(df1, 0,1,200,Y0,ud1,x);
    // //df1 funkcja ktora zwraca rownanie rozniczkowe
    // int n=get_len(Y[0]);
    // double max=Y[1](0,2); //(0,2)? -> tak jest skonstruowana matrix
    // for(int i=0;i<n;i++){
    //     if(max<Y[1](i,2))
    //         max=Y[1](i,2);
    // }
    // y=abs(max-50);
    // return y;
}