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
    //dla tego nie uzywamy metody ekspansjif1R
    matrix y;
    matrix Y0=matrix(3,new double[]{5,1,20});
    matrix* Y= solve_ode(df1, 0,1,2000,Y0,ud1,x);
    //df1 funkcja ktora zwraca rownanie rozniczkowe
    int n=get_len(Y[0]);
    double max=Y[1](0,2); //(0,2)? -> tak jest skonstruowana matrix
    for(int i=0;i<n;i++){
        if(max<Y[1](i,2))
            max=Y[1](i,2);
    }
    y=abs(max-50);

    return y;
}

matrix df1(double t, matrix Y, matrix ud1, matrix ud2) {
	// wektor zmian po czasie
	matrix dY(3,1);

	// zmienne zadania
	double Va = Y(0);
	double Vb = Y(1);
	double Tb = Y(2);
	// dane z zadania
	double Pa = ud1(0);
	double Ta = ud1(1);
	double Pb = ud1(2);
	double Db = ud1(3);
	double F_in = ud1(4);
	double T_in = ud1(5);
	double a = ud1(6);
	double b = ud1(7);
	double g = ud1(8);
	double Da = ud2(0);	// to bedziemy optymalizowac
	// std::cout << ud2 << endl;
	// obliczanie wylanego strumienia ze zbiornika A

	double Fa_out{};
	if (Va > 0.0)
		Fa_out = a * b * Da * sqrt(2 * g * Va / Pa); // Va to nasz Y(0)
	// obliczanie wylanego strumienia ze zbiornika B
	double Fb_out{};
	if (Vb > 0.0)
		Fb_out = a * b * Db * sqrt(2 * g * Vb / Pb); // Vb to nasz Y(1)
	// double dTb_dt = Fa_out/Vb * (Ta - Tb) + F_in/Vb * (T_in - Tb);	// Tb to nasz Y(2)
	//Ustalanie zmiany obj�to�ci w zbiornku A
	if (Y(0) + dY(0) < 0)
		dY(0) = -Y(0); //Wylanie reszty je�li strumie� wi�kszy od obj�to�ci wody
	else
		dY(0) = -Fa_out; //Wylanie strumienia
	//Ustalanie zmien obj�to�ci w zbiorniku B
	if (Y(1) + dY(1) < 0)
		dY(1) = -Y(1); //Wylanie reszty je�li strumie� wi�kszy od obj�to�ci wody
	else
		dY(1) = Fa_out + F_in - Fb_out; //Wylanie strumienia oraz wlanie wody z kranu i ze zbiornika A
	//Ustalenie zmian temperatury w zbiorniku B
	if (Vb > 0)
		dY(2) = (F_in / Vb) * (T_in - Tb) + (Fa_out / Vb) * (Ta - Tb); //Formu�a je�li zbiornik B nie jest pusty
	else
		dY(2) = 0; //Pusty zbiornik B
	//Zwracanie zmian po czasie
	return dY;
}

// lab 2 functions
matrix ff2T(matrix x1, matrix x2, matrix ud1, matrix ud2) {
	matrix y;
	y = pow(m2d(x1), 2) + pow(m2d(x2), 2) - cos(2.5 * M_PI * m2d(x1)) - cos(2.5 * M_PI * m2d(x2)) + 2;
	return y;
}

matrix ff2R(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	y = 0;

	// warunki poczatkowe
	matrix Y0(2,1);
	matrix Yref(2, new double[]{3.14, 0});

	Y0 = ud1(4);
	Y0 = ud1(5);

	Yref = Y0(1);
	// Yref =

	matrix * Y = solve_ode(df, 0, 0.1, 100, Y0, Yref, x);

	double alpha_ref = ud1(4);
	double omega_ref = ud1(5);

	y = 0;
	int n = get_len(Y[0]);

	for(int i = 0; i < n; i++) {
	// wedlug funkcji podcalkowej
		y=y+10*pow(Yref(0) - Y[1](i,0),2) + pow*(Yref(1) - Y[1](i,1),2) + pow(x(0) *( Yref(0)-Y[1](i,0)) +
			+ x(1) *(Yref(1)-Y[1](i,1) ),2);
	}

	y = y * 0.1;


}