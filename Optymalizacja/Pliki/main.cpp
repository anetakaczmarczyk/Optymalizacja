/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include"opt_alg.h"

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		lab1();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	system("pause");
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;
	int Nmax = 10000;
	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);  //lb dolne ograniczenie, ub gorne ograniczenie, od -5 do 5
	solution opt;
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Wahadlo
	Nmax = 1000;
	epsilon = 1e-2;
	lb = 0;
	ub = 5;
	double teta_opt = 1;
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(opt.x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
	ofstream Sout("symulacja_lab0.csv");
	Sout << hcat(Y[0], Y[1]);
	Sout.close();
	Y[0].~matrix();
	Y[1].~matrix();
}

void lab1()
{
	// solution do testow
	solution test_opt;

	//martynka expension test
	double epsilon = 1e-18;
	double d = 0.01;
	int Nmax = 200;
	double alpha = 1.1;


	// Generator losowania liczb
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> x0_dist(0.0, 100.0);

	std::stringstream test_ss;

	double* result=nullptr;

	for (int i = 0; i < 100; ++i)
	{
		double x0 = x0_dist(gen);
		result = expansion(ff1, x0, d, alpha, Nmax);

		// Zapis wyników do stringstream
		test_ss << x0 << ";" << result[0] << ";" << result[1] << ";" << result[2] << ";\n";
		solution::clear_calls();

		// obliczanie minimum metoda fibonacciego
		test_opt = fib(ff1, result[0], result[1], epsilon);
		// zapis do stringa: x_min; y_min; f_calls;
		test_ss << m2d(test_opt.x) << "; " << m2d(test_opt.y) << "; " << test_opt.f_calls << ";\n";
		solution::clear_calls();

		// Zwolnienie pamięci
		delete[] result;
	}
	std::cout << "Wyniki testu ekspansji:\n";
	std::cout << test_ss.str() << std::endl;

	solution::clear_calls();

	// ania - fibonacci w przedziale [-100; 100]
	test_opt = fib(ff1, -100 , 100 , epsilon);
	std::cout << "Minimum metoda Fibonacci'ego:\n";
	std::cout << test_opt << std::endl;
	solution::clear_calls();
}

void lab2()
{
	//Aneta test
	int x = 2;
}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
