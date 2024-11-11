#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
    try
    {
        int i = 0;
        solution xi_sol, xi_next_sol;
        double xi, xi_next;

        xi = x0;
        xi_next = xi + d;

        xi_sol.x = xi;
        xi_sol.fit_fun(ff, ud1, ud2);

        xi_next_sol.x = xi_next;
        xi_next_sol.fit_fun(ff, ud1, ud2);


        if (xi_next_sol.y == xi_sol.y)
            return new double[3]{ xi, xi_next, (double)solution::f_calls };

        if (xi_next_sol.y > xi_sol.y)
        {
            d = -d;
            xi_next = xi + d;
            xi_next_sol.x = xi_next;
            xi_next_sol.fit_fun(ff, ud1, ud2);


            if (xi_next_sol.y >= xi_sol.y)
                return new double[3]{ xi_next, xi - d, (double)solution::f_calls };
        }

        solution::clear_calls();
        double xi_prev{};
        double f_xi = m2d(xi_sol.y);
        do
        {
            if (solution::f_calls > Nmax)
            {
                xi_next_sol.flag = 0;
                throw std::string("Maximum amount of f_calls reached!");
            }

            ++i;
            xi_next = xi + pow(alpha, i) * d;

            xi_next_sol.x = xi_next;
            xi_next_sol.fit_fun(ff, ud1, ud2);

            if (!(f_xi > xi_next_sol.y))
                break;

            xi_prev = xi;
            xi = xi_next;
            f_xi = m2d(xi_next_sol.y);

        } while (true);

        if (d > 0)
            return new double[3]{ xi_prev, xi_next, (double)solution::f_calls };

        return new double[3]{ xi_next, xi_prev, (double)solution::f_calls };
    }
    catch (const std::string& ex_info)
    {
        throw ("double* expansion(...):\n" + ex_info);
    }
}


solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		// ustalenie końców początkowego przedziału poszukiwań
		std::vector<double> phi;
		phi.push_back(0);
		phi.push_back(1);

		// poszukiwanie najmniejszego k, dla ktorego phi_k > L/epsilon
		int k = 1;
		while (phi[k] <= (b - a) / epsilon){
			phi.push_back(phi[k] + phi[k - 1]);
			k++;
		}

		double a0 = a;
		double b0 = b;
		double c0 = b0 - phi[k-1] / phi[k] * (b0 - a0);
		double d0 = a0 + b0 - c0;

		std::stringstream fib_ss;	// do zapisu danych

		solution c_sol, d_sol;	// jak sie chce miec wartosc funkcji w punkcie to trzeba uzyc klasy solution xd
		for (int i = 0; i <= k - 3; ++i)
		{
			fib_ss << b0-a0 << ";" << std::endl;
			c_sol.x = c0;
			c_sol.fit_fun(ff, ud1);	// to uzupelnia pole y w klasie solution

			d_sol.x = d0;
			d_sol.fit_fun(ff, ud1);

			// redukowanie przedziału
			if (c_sol.y < d_sol.y) {
				b0 = d0;
			}
			else {
				a0 = c0;
			}

			c0 = b0 - phi[k - i - 2] / phi[k - i - 1] * (b0 - a0);
			d0 = a0 + b0 - c0;
		}

		// zapis wynikow do pliku
		// std::ofstream file(R"(C:\Users\Ania\CLionProjects\Optymalizacja\Optymalizacja\lab1-analiza\lab1-100-przedzialow-fib.txt)"); //musialam dac cala sciezke bo nie dzialalo xd
		// if (file.is_open()) {
		// 	file << fib_ss.str();
		// 	file.close();
		// }else {
		// 	cerr << "Nie udało się otworzyć pliku do zapisu.\n";
		// }

		// zwracanie rozwiazania
		solution Xopt;
		Xopt.x = c0;
		Xopt.fit_fun(ff, ud1);

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		int i = 0;
		double d_prev{};
		double d{};
		double c = (a+b)/2;
		solution a_sol, b_sol, c_sol, d_sol;
		std::stringstream lag_ss;	// do zapisu danych
		do {
			lag_ss << b-a << ";" << std::endl;
			d_prev = d;
			a_sol.x = a;
			b_sol.x = b;
			c_sol.x = c;
			a_sol.fit_fun(ff, ud1);
			b_sol.fit_fun(ff, ud1);
			c_sol.fit_fun(ff, ud1);

			double l = m2d(a_sol.y) * (pow(b,2) - pow(c, 2));
			l += m2d(b_sol.y) * (pow(c, 2) - pow(a, 2));
			l += m2d(c_sol.y) * (pow(a, 2) - pow(b, 2));
			double m = m2d(a_sol.y) * (b - c);
			m += m2d(b_sol.y) * (c - a);
			m += m2d(c_sol.y) * (a - b);

			if (m<=0) {
				Xopt.flag = 0;
				break;
			}
			d = 0.5 * (l / m);
			d_sol.x = d;
			d_sol.fit_fun(ff, ud1);
			if (a<d && d < c) {
				if (d_sol.y < c_sol.y) {
					//a = a;
					b = c;
					c = d;
				}
				else {
					a = d;
					//c = c;
					//b = b;
				}
			}
			else if (c < d && d < b) {
				if (d_sol.y < c_sol.y) {
					a = c;
					c = d;
					//b = b;
				}
				else {
					//a = a;
					//c = c;
					b = d;
				}
			}
			else {
				Xopt.flag = 0;
				break;
			}
			i+= 1;
			if (solution::f_calls > Nmax) {
				Xopt.flag = 0;
				break;
			}

		}
		while (!(b - a < epsilon or abs(d - d_prev) <= gamma));
		// zapis wynikow do pliku
		// std::ofstream file(R"(C:\Users\Animatt\CLionProjects\Optymalizacja\Optymalizacja\lab1-analiza\lab1-100-przedzialow-lag.txt)"); //musialam dac cala sciezke bo nie dzialalo xd
		// if (file.is_open()) {
		// 	file << lag_ss.str();
		// 	file.close();
		// }else {
		// 	cerr << "Nie udało się otworzyć pliku do zapisu.\n";
		// }
		Xopt.x = d;
		Xopt.fit_fun(ff, ud1);
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{//etap probny -> badamy lokalne zachowanie funkcji wokol[ [unktu bazowego i uzywamy tego punktu bazowaego
    //etap roboczy -> przenosimys sie do nowego punktu i w nim robimy znowu etap probny

	try
	{

		solution Xopt;
		//Tu wpisz kod funkcji
		solution XB, X;
		X.x = x0;
		while (s<epsilon) {
			XB.x = X.x;
			XB.fit_fun(ff, ud1);
			X = HJ_trial(ff, XB, s, ud1, ud2);
			X.fit_fun(ff, ud1);
			if(X.y < XB.y) {
				while(X.y < XB.y) {
					solution tempXB = XB;
					XB = X;
					X.x = 2*XB.x - tempXB.x;
					X = HJ_trial(ff, X, s, ud1, ud2);
					X.fit_fun(ff, ud1);
					if (solution::f_calls > Nmax) {
						break;
					}
				}
			}
			else {
				s = alpha * s;
			}
			if(solution::f_calls > Nmax) {
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		int n = 2; // Bo 2 wymiary, byłoby 3 jakby 3 wymiary

		//Tu wpisz kod funkcji
		for (int i = 0; i<n; i++) {
			XB.fit_fun(ff, ud1);
			solution f1_sol, f2_sol;
			f1_sol.x = XB.x + s;
			f1_sol.fit_fun(ff, ud1);
			f2_sol.x = XB.x - s;
			f2_sol.fit_fun(ff, ud1);
			if (f1_sol.y < XB.y) {
				XB.x = m2d(XB.x) + s;
			}
			else if(f2_sol.y < XB.y) {
				XB.x = m2d(XB.x) - s;
			}
		}
		return XB.x;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{//wektor p to wektor porazek jak moje zycie
    //i nr czegos i J nr bazu? nwm do ktorego alg to bo nie sluchalam
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji
        //macierz D to macierz kierunkow  d1 d2 i potem macierz
        //zlozona z lambd
        //zaimplementowac w 2D


		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
