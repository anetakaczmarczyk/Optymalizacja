/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/
#include <iostream>
#include <filesystem> // C++17 i nowsze
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
		lab2();
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

	double epsilon = 1e-18;
	double gamma = 1e-30;
	double d = 0.01;
	int Nmax = 200;
	double alpha = 1.1;

	std::stringstream test_ss;	// do zapisu danych

	// ekspansja dla 3 roznych wspolczynnikow
	for (int j = 0; j < 3; j++) {
		// Generator losowania liczb
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> x0_dist(0.0, 100.0);

		double* result=nullptr;

		for (int i = 0; i < 100; ++i)
		{
			double x0 = x0_dist(gen);
			result = expansion(ff1, x0, d, alpha, Nmax);
			// Zapis wyników do stringstream
			test_ss << x0 << ";" << result[0] << ";" << result[1] << ";" << result[2] << ";";
			solution::clear_calls();

			// obliczanie minimum metoda fibonacciego
			test_opt = fib(ff1, result[0], result[1], epsilon);
			// zapis do stringa: x_min; y_min; f_calls; lokalne/globalne;
			test_ss << m2d(test_opt.x) << ";" << m2d(test_opt.y) << ";" << test_opt.f_calls << ";" << (test_opt.x > -1 && test_opt.x < 1 ? "lokalne" : "globalne") << ";";
			solution::clear_calls();

			// obliczanie minimum metoda lagrange'a
			test_opt = lag(ff1, result[0], result[1], epsilon, gamma, Nmax);
			// zapis do stringa: x_min; y_min; f_calls
			test_ss << m2d(test_opt.x) << ";" << m2d(test_opt.y) << ";" << test_opt.f_calls << ";" << (test_opt.x > -1 && test_opt.x < 1 ? "lokalne" : "globalne") << ";\n";

			// Zwolnienie pamięci
			delete[] result;
		}
		alpha *= 2;
	}

	// zapis wynikow do pliku
	// std::ofstream file("C:\\Users\\Ania\\CLionProjects\\Optymalizacja\\Optymalizacja\\lab1-analiza\\lab1-100-wynikow.txt"); //musialam dac cala sciezke bo nie dzialalo xd
	// if (file.is_open()) {
	// 	file << test_ss.str();
	// 	file.close();
	// }else {
	// 	cerr << "Nie udało się otworzyć pliku do zapisu.\n";
	// }

	std::cout << "Wyniki:\n";
	std::cout << test_ss.str() << std::endl;

	// przyklad bez eskpansji [-100; 100]
	test_opt = fib(ff1, -100 , 100 , epsilon);
	std::cout << "Minimum metoda Fibonacci'ego:\n";
	std::cout << test_opt << std::endl;
	solution::clear_calls();

	test_opt = lag(ff1, -100, 100, epsilon, gamma, Nmax);
	std::cout << "Minimum metoda Lagrange'a:\n";
	std::cout << test_opt << std::endl;
	solution::clear_calls();

	// problem rzeczywisty

	matrix ud1 = matrix(9, 1);
	ud1(0, 0) = 0.5;	// Pa - pole podst. A
	ud1(1, 0) = 90;		// Ta - temperatura
	ud1(2, 0) = 1;		// Pb - pola podst. B
	ud1(3, 0) = 0.00365665;	// Db (m^2) - pole przekroju
	ud1(4, 0) = 0.01;	//FIN (m^3/s) - szybkosc wlewania
	ud1(5, 0) = 20;		// T_in - temperatura
	ud1(6, 0) = 0.98;	// a - lepkosc
	ud1(7, 0) = 0.63;	// b - zwezenie strumienia
	ud1(8, 0) = 9.81;	// g - przysp. graw.
	//Zakres szukania Da
	double Da_0_s = 1.0 * 0.0001;
	double Da_0_f = 100 * 0.0001;
	solution::clear_calls();

	//Szukanie minimum metoda fibonacciego
	solution opt = fib(f1R, Da_0_s, Da_0_f, epsilon, ud1);
	std::cout << "minimum fib dla rzeczywistego\n" << opt;
	solution::clear_calls();

	//Szukanie minimum metoda lagrange'a
	opt = lag(f1R, Da_0_s, Da_0_f, epsilon, gamma, Nmax, ud1);
	std::cout << "minimum lag dla rzeczywistego\n" << opt;
	solution::clear_calls();

	//Warunki pocz¹tkowe
	matrix Y0 = matrix(3, 1);
	Y0(0) = 5.0; //Poczatkowa objetosc w a
	Y0(1) = 1.0; //Poczatkowa objetosc w b
	Y0(2) = 20.0;//Poczatkowa temperatura w b

	//Symulacja
	matrix* Y = solve_ode(df1, 0, 1, 2000, Y0, ud1, opt.x);
	std::stringstream symulationFib_ss;	// do zapisu danych
	symulationFib_ss << hcat(Y[0], Y[1]) << ";";
	// zapis wynikow do pliku
	std::ofstream file(R"(C:\Users\Animatt\CLionProjects\Optymalizacja\Optymalizacja\lab1-analiza\lab1-symulation-fib.txt)"); //musialam dac cala sciezke bo nie dzialalo xd
	if (file.is_open()) {
		file << symulationFib_ss.str();
		file.close();
	}else {
		cerr << "Nie udało się otworzyć pliku do zapisu.\n";
	}


	// Max temp
	int n = get_len(Y[0]);
	double Tb_max = Y[1](0, 2);
	for (int i = 0; i < n; ++i)
	{
		if (Tb_max < Y[1](i, 2))
			Tb_max = Y[1](i, 2);
	}
	std::cout << "Fibonacci tb_max: " << Tb_max << "\n\n";

	//Szukanie minimum
	opt = lag(f1R, Da_0_s, Da_0_f, epsilon, epsilon, Nmax, ud1);
	std::cout << opt;
	solution::clear_calls();

	//Symulacja
	Y = solve_ode(df1, 0, 1, 2000, Y0, ud1, opt.x);
	std::stringstream symulationLag_ss;	// do zapisu danych
	symulationLag_ss << hcat(Y[0], Y[1]) << ";";
	// zapis wynikow do pliku
	std::ofstream file1(R"(C:\Users\Ania\CLionProjects\Optymalizacja\Optymalizacja\lab1-analiza\lab1-symulation-lagtxt)"); //musialam dac cala sciezke bo nie dzialalo xd
	if (file1.is_open()) {
		file1 << symulationLag_ss.str();
		file1.close();
	}else {
		cerr << "Nie udało się otworzyć pliku do zapisu.\n";
	}
	// Max temp
	 n = get_len(Y[0]);
	 Tb_max = Y[1](0, 2);
	for (int i = 0; i < n; ++i)
	{
		if (Tb_max < Y[1](i, 2))
			Tb_max = Y[1](i, 2);
	}
	std::cout << "Lagrange tb_max: " << Tb_max << "\n\n";
}

void lab2()
{
	//warunki do obu
	double s = 0.1;
	double alpha = 0.2;
	double beta = 0.2;
	double epsilon = 1E-6;
	int Nmax = 2000;
	double alphaRosen = 1.2;

	double tolerance = 0.01;
	solution test_opt;
	std::stringstream test_ss;	// do zapisu danych

	// dla 3 różnych kroków s: 0.1, 1.1, 2.1
	for (int j = 0; j < 3; j++) {
		// Generator losowania liczb
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> x0_dist(-1.0, 1.0);

		for (int i = 0; i < 100; ++i)
		{
			matrix x0 = matrix(2, new double[2] {x0_dist(gen), x0_dist(gen)});
			//zapis do stringa: x1; x2 wygenerowane
			test_ss << x0(0) << ";" << x0(1) << ";";
			test_opt = HJ(ff2T, x0,  s, alpha, epsilon, Nmax);

			// zapis do stringa: x1_min; x2_min; y_min; f_calls; globalne/lokalne;
			test_ss << m2d(test_opt.x(0)) << ";"<< m2d(test_opt.x(1)) << ";"<< m2d(test_opt.y) << ";" << test_opt.f_calls << ";" << (abs(m2d(test_opt.y) ) < tolerance ? "TAK" : "NIE") << ";";
			solution::clear_calls();

			test_opt = Rosen(ff2T, x0, matrix(2, new double[2] {s, s}), alphaRosen, beta, epsilon, Nmax);
			// zapis do stringa: x1_min; x1_min; y_min; f_calls; globalne/lokalne
			test_ss << m2d(test_opt.x(0)) << ";"<< m2d(test_opt.x(1)) << ";"<< m2d(test_opt.y) << ";" << test_opt.f_calls << ";" << (abs(m2d(test_opt.y) ) < tolerance ? "TAK" : "NIE") << ";\n";
			solution::clear_calls();
		}
		s*=2;
	}

	// zapis wynikow do pliku
	std::ofstream file("Optymalizacja\\lab2-analiza\\lab2-100-optymalizacji.txt");
	if (file.is_open()) {
		file << test_ss.str();
		file.close();
	}else {
		cerr << "Nie udało się otworzyć pliku do zapisu.\n";
	}

	std::cout << "Wyniki:\n";
	std::cout << test_ss.str() << std::endl;

	//powrót do kroku równego 0.1
	s = 0.1;
	//Excel - wykres
	// matrix x0 = matrix(2, new double[2] {-0.45, 0.45});
	// HJ(ff2T, x0,  s, alpha, epsilon, Nmax);
	// Rosen(ff2T, x0, matrix(2, new double[2] {s, s}), alphaRosen, beta, epsilon, Nmax);

  
	// Problem rzeczywisty
	// dane do zadania
	matrix ud1(6, 1);
	ud1(0, 0) = 1.0;	// dlugosc ramienia (m)
	ud1(1, 0) = 1.0;	// masa ramienia (kg)
	ud1(2, 0) = 5.0;	// masa ciezarku (kg)
	ud1(3, 0) = 0.5;	// wspolczynnik tarcia (Nms)
	ud1(4, 0) = 3.14;	// alpha_ref (rad)
	ud1(5, 0) = 0.0;	// omega_ref (rad/s)

	// warunki poczatkowe
	matrix k_0 (2, 1);
	k_0(0, 0) = 1.0;
	k_0(1, 0) = 1.0;

	matrix Y = matrix(2, 1);

	solution realHJ = HJ(ff2R, k_0, s, alpha, epsilon, Nmax,  ud1);
	std::cout << "HJ realistyczny\n" << realHJ;
	solution::clear_calls();
	// matrix* X = solve_ode(df2, 0.0, 0.1, 100.0, Y, ud1, realHJ.x);

	// std::stringstream symulationHJ_ss;	// do zapisu danych
	// symulationHJ_ss << hcat(X[0], X[1]) << ";";
	// // zapis wynikow do pliku
	// std::ofstream file(R"(C:\Users\Animatt\CLionProjects\Optymalizacja\Optymalizacja\lab2-analiza\lab2-symulation-HJ.txt)"); //musialam dac cala sciezke bo nie dzialalo xd
	// if (file.is_open()) {
	// 	file << symulationHJ_ss.str();
	// 	file.close();
	// }else {
	// 	cerr << "Nie udało się otworzyć pliku do zapisu.\n";
	// }

	// solution realRosen = Rosen(ff2R, k_0, matrix(2, new double[2] {s, s}), alphaRosen, beta, epsilon, Nmax, ud1);
	// std::cout << "Rosen realistyczny\n" << realRosen;
	// solution::clear_calls();
	// matrix* X = solve_ode(df2, 0.0, 0.1, 100.0, Y, ud1, realRosen.x);
	// std::stringstream symulationRosen_ss;	// do zapisu danych
	// symulationRosen_ss << hcat(X[0], X[1]) << ";";
	// // zapis wynikow do pliku
	// std::ofstream file(R"(C:\Users\Animatt\CLionProjects\Optymalizacja\Optymalizacja\lab2-analiza\lab2-symulation-Rosen.txt)"); //musialam dac cala sciezke bo nie dzialalo xd
	// if (file.is_open()) {
	// 	file << symulationRosen_ss.str();
	// 	file.close();
	// }else {
	// 	cerr << "Nie udało się otworzyć pliku do zapisu.\n";
	// }
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
