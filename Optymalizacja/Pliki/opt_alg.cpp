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

		solution c_sol, d_sol;	// jak sie chce miec wartosc funkcji w punkcie to trzeba uzyc klasy solution xd
		for (int i = 0; i <= k - 3; ++i)
		{
			c_sol.x = c0;
			c_sol.fit_fun(ff);	// to uzupelnia pole y w klasie solution

			d_sol.x = d0;
			d_sol.fit_fun(ff);

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

		// zwracanie rozwiazania
		solution Xopt;
		Xopt.x = c0;
		Xopt.fit_fun(ff);

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
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

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
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

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
