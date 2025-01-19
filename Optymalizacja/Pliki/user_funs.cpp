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
matrix ff2T(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2;
	return y;
}
matrix df2(double t, matrix Y, matrix ud1, matrix ud2) {
	matrix dY(2, 1);
	double alpha = Y(0);
	double omega = Y(1);

	double l = ud1(0);
	double m_r = ud1(1);
	double m_c = ud1(2);
	double b = ud1(3);
	double alpha_ref = ud1(4);
	double omega_ref = ud1(5);

	double k1 = ud2(0);
	double k2 = ud2(1);

	double I = (1.0 / 3.0) * m_r * pow(l, 2) + m_c * pow(l, 2);

	double Mt = k1 * (alpha_ref - alpha) + k2 * (omega_ref - omega);
	dY(0) = Y(1);
	dY(1) = (Mt -b * omega)/ I;
	return dY;
}
matrix ff2R(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	y = 0;

	// warunki poczatkowe
	matrix Y0(2,1);

	matrix * Y = solve_ode(df2, 0, 0.1, 100, Y0, ud1, x);

	double alpha_ref = ud1(4);
	double omega_ref = ud1(5);

	int n = get_len(Y[0]);

	for(int i = 0; i < n; i++) {
	// wedlug funkcji podcalkowej
		y=y+10*pow(alpha_ref - Y[1](i,0),2) + pow(omega_ref - Y[1](i,1),2) + pow(x(0) *( alpha_ref-Y[1](i,0)) +
			+ x(1) *(omega_ref-Y[1](i,1) ),2);
	}

	y = y * 0.1;

	return y;
}

// lab 3 functions
matrix ff3T_out(matrix x, matrix ud1, matrix ud2) {
	matrix y;;
	y = sin(M_PI * sqrt(pow(x(0)/M_PI, 2) + pow(x(1)/M_PI, 2)))/M_PI * sqrt(pow(x(0)/M_PI, 2) + pow(x(1)/M_PI, 2));
	if( -x(0) + 1 > 0) {
		y = y + ud2 * pow(-x(0) + 1, 2);
	}
	if( -x(1) + 1 > 0) {
		y = y + ud2 * pow(-x(1) + 1, 2);
	}
	if( norm(x) - ud1 > 0) {
		y = y + ud2 * pow(norm(x) - ud1, 2);
	}
	return y;
}
matrix ff3T_in(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	y = sin(M_PI * sqrt(pow(x(0)/M_PI, 2) + pow(x(1)/M_PI, 2)))/M_PI * sqrt(pow(x(0)/M_PI, 2) + pow(x(1)/M_PI, 2));
	if( -x(0)+1 >=0) {
		y=1e10;
	}else {
		y=y-ud2/(-x(0)+1);
	}
	if( -x(1)+1 >=0) {
		y=1e10;
	}else {
		y=y-ud2/(-x(1)+1);
	}
	if( norm(x) - ud1 >0) {
		y=1e10;
	}else {
		y=y-ud2/(norm(x) - ud1);
	}
	return y;
}

matrix df3(double t, matrix Y, matrix ud1, matrix ud2) {
	matrix dY(4,1);
	double C = ud1(0);
	double rho = ud1(1);
	double r = ud1(2);
	double m = ud1(3);
	double g = ud1(4);

	double s = M_PI * pow(r, 2);
	double Dx = 0.5 * C * rho * s * Y(1) * abs(Y(1));
	double Dy = 0.5 * C * rho * s * Y(3) * abs(Y(3));
	double Fmx = rho * Y(3) * ud2(0) * M_PI * pow(r, 3);
	double Fmy = rho * Y(1) * ud2(0) * M_PI * pow(r, 3);

	dY(0) = Y(1);
	dY(1) = (-Dx - Fmx) / m;
	dY(2) = Y(3);
	dY(3) = ((-m * g) - Dy - Fmy) / m;

	return dY;
}

matrix ff3R(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	matrix Y0(4, new double[4] {0.0, x(0), 100.0, 0.0});
	matrix* Y = solve_ode(df3, 0.0, 0.01, 7.0, Y0, ud1, x(1));
	int n = get_len(Y[0]);
	int i0 = 0;
	int i50 = 0;

	for (int i = 0; i < n; i++)
	{
		if (abs(Y[1](i, 2) - 50.0) < abs(Y[1](i50, 2) - 50.0))
			i50 = i;
		if (abs(Y[1](i, 2)) < abs(Y[1](i0, 2)))
			i0 = i;
	}

	y = -Y[1](i0, 0);

	if (abs(x(0)) - 10 > 0)
		y = y + ud2 * pow(abs(x(0)) - 10, 2);
	if (abs(x(1)) - 15 > 0)
		y = y + ud2 * pow(abs(x(1)) - 15, 2);
	if (abs(Y[1](i50, 0) - 5.0) - 0.5 > 0)
		y = y + ud2 * pow(abs(Y[1](i50, 0) - 5.0) - 0.5, 2);

	return y;
}

matrix ff4T(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	if (isnan(ud2(0, 0))) {
		y = pow(x(0) + 2 * x(1) - 7, 2) + pow(2*x(0)+x(1)-5,2);
	}else {
		y = ff4T(ud2[0] + x * ud2[1]);
	}
	return y;
}

matrix gfT(matrix x, matrix ud1, matrix ud2) {
	matrix y(2, 1);
	y(0) = -34.0 + 10.0 * x(0) + 8.0 * x(1);
	y(1) = -34.0 + 8.0 * x(0) + 10.0 * x(1);
	return y;
}

matrix hfT(matrix x, matrix ud1, matrix ud2) {
	matrix y(2, 2);
	y(0, 0) = 10;
	y(0, 1) = 8;
	y(1, 0) = 8;
	y(1, 1) = 10;
	return y;
}

double sygmoid(matrix theta, matrix x) {
	return 1.0/(1.0 + exp(-1.0 * m2d(trans(theta)*x)));
}

matrix ff4R(matrix theta, matrix X, matrix Y) {
	matrix y;
	double sum = 0.0;
	for (int i = 0; i < 100; i++) {
		double yi = Y[i](0);
		matrix xi = X[i];
		sum += yi * log(sygmoid(theta, xi)) + (1.0 - yi) * log(1.0 - sygmoid(theta, xi));
	}
	y = (-1.0/100.0) * sum;
	return y;
}

matrix gf4R(matrix theta, matrix X, matrix Y) {
	matrix y(3,1);
	for(int i = 0; i<3; i++) {
		double sum = 0.0;
		for (int j =0; j<100; j++) {
			double yi = Y[j](0);
			matrix xi = X[j];
			sum += (sygmoid(theta, xi) - yi) * xi(i);
		}
		y(i) = (1.0/100.0) * sum;
	}
	return y;
}

matrix ff5T(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	if (isnan(ud2(0, 0))) {
		y=matrix(2, 1);
		y(0) = ud1(1)* (pow(x(0)-2, 2) + pow(x(1) - 2, 2));
		y(1) = 1/ud1(1) * (pow(x(0)+2, 2) + pow(x(1) + 2, 2));

	} else {
		matrix yt;
		yt = ff5T(ud2[0]+x*ud2[1], ud1);
		y = ud1(0) * yt(0) + (1-ud1(0))*yt(1);
	}
	return y;
}

matrix ff5R(matrix x, matrix ud1, matrix ud2) {
	matrix y;

	double rho = 7800;
	double P = 1000;
	double E = 207e9;

	if (isnan(ud2(0, 0))) {
		double l = x(0);
		double d = x(1);
		y = matrix(3, 1);
		y(0) = rho * l * M_PI * (pow(d,2) / 4);
		y(1) = (64 * P * pow(l, 3)) / (3 * E * M_PI * pow(d, 4));
		y(2) = (32 * P * l)/(M_PI * pow(d, 3));
	}
	else
	{
		matrix xt = ud2[0] + x * ud2[1];

		matrix yt = ff5R(xt, ud1);
		y = ud1 * (yt(0) - 0.12) / (15.3 - 0.12) + (1 - ud1) * (yt(1) - 4.2e-5) / (3.28 - 4.2e-5);

		const double c = 1e10;

		if (xt(0) < 0.2) y = y + c * pow(0.2 - xt(0), 2);
		if (xt(0) > 1) y = y + c * pow(xt(0) - 1, 2);
		if (xt(1) < 0.01) y = y + c * pow(0.01 - xt(1), 2);
		if (xt(1) > 0.05) y = y + c * pow(xt(1) - 0.05, 2);
		if (xt(1) > 0.005) y = y + c * pow(xt(1) - 0.005, 2);

		if (yt(1) > 0.005) y = y + c * pow(yt(1) - 0.005, 2);
		if (yt(2) > 300e6) y = y + c * pow(yt(2) - 300e6, 2);
	}
	return y;
}

matrix ff6T(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2;
	return y;
}

matrix df6(double t, matrix Y, matrix ud1, matrix ud2) {
	matrix dY(4,1);
	double m1 = 5;
	double m2 = 5;
	double k1 = 1;
	double k2 = 1;

	double b1 = ud2(0);
	double b2 = ud2(1);

	double F = 1;

	dY(0) = Y(1);
	dY(1) = (-b1 * Y(1) - b2 * (Y(1) - Y(3)) - k1 * Y(0) - k2 * (Y(0) - Y(2))) / m1;
	dY(2) = Y(3);
	dY(3) = (F + b2 * (Y(1) - Y(3)) + k2 * (Y(0) - Y(2))) / m2;

	return dY;
}

matrix ff6R(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	matrix Y0(4, 1);
	matrix* Y = solve_ode(df6, 0, 0.1, 100, Y0, ud1, x[0]);
	for (int i = 0; i < ud1(0); i++)
	{
		y = y + abs(ud2(i, 0) - Y[1](i, 0)) + abs(ud2(i, 1) - Y[1](i, 2));
	}
	y(0) = y(0) / (2 * ud1(0));


	return y;
}



