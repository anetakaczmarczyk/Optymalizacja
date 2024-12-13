#include"opt_alg.h"

solution MC(matrix (*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1,
            matrix ud2) {
    try {
        solution Xopt;
        while (true) {
            Xopt = rand_mat(N);
            for (int i = 0; i < N; ++i)
                Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
            Xopt.fit_fun(ff, ud1, ud2);
            if (Xopt.y < epsilon) {
                Xopt.flag = 1;
                break;
            }
            if (solution::f_calls > Nmax) {
                Xopt.flag = 0;
                break;
            }
        }
        return Xopt;
    } catch (string ex_info) {
        throw ("solution MC(...):\n" + ex_info);
    }
}

double *expansion(matrix (*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1,
                  matrix ud2) {
    try {
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
            return new double[3]{xi, xi_next, (double) solution::f_calls};

        if (xi_next_sol.y > xi_sol.y) {
            d = -d;
            xi_next = xi + d;
            xi_next_sol.x = xi_next;
            xi_next_sol.fit_fun(ff, ud1, ud2);


            if (xi_next_sol.y >= xi_sol.y)
                return new double[3]{xi_next, xi - d, (double) solution::f_calls};
        }

        solution::clear_calls();
        double xi_prev{};
        double f_xi = m2d(xi_sol.y);
        do {
            if (solution::f_calls > Nmax) {
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
            return new double[3]{xi_prev, xi_next, (double) solution::f_calls};

        return new double[3]{xi_next, xi_prev, (double) solution::f_calls};
    } catch (const std::string &ex_info) {
        throw ("double* expansion(...):\n" + ex_info);
    }
}


solution fib(matrix (*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2) {
    try {
        // ustalenie końców początkowego przedziału poszukiwań
        std::vector<double> phi;
        phi.push_back(0);
        phi.push_back(1);

        // poszukiwanie najmniejszego k, dla ktorego phi_k > L/epsilon
        int k = 1;
        while (phi[k] <= (b - a) / epsilon) {
            phi.push_back(phi[k] + phi[k - 1]);
            k++;
        }

        double a0 = a;
        double b0 = b;
        double c0 = b0 - phi[k - 1] / phi[k] * (b0 - a0);
        double d0 = a0 + b0 - c0;

        std::stringstream fib_ss; // do zapisu danych

        solution c_sol, d_sol; // jak sie chce miec wartosc funkcji w punkcie to trzeba uzyc klasy solution xd
        for (int i = 0; i <= k - 3; ++i) {
            fib_ss << b0 - a0 << ";" << std::endl;
            c_sol.x = c0;
            c_sol.fit_fun(ff, ud1); // to uzupelnia pole y w klasie solution

            d_sol.x = d0;
            d_sol.fit_fun(ff, ud1);

            // redukowanie przedziału
            if (c_sol.y < d_sol.y) {
                b0 = d0;
            } else {
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
    } catch (string ex_info) {
        throw ("solution fib(...):\n" + ex_info);
    }
}

solution lag(matrix (*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax,
             matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        int i = 0;
        double d_prev{};
        double d{};
        double c = (a + b) / 2;
        solution a_sol, b_sol, c_sol, d_sol;
        std::stringstream lag_ss; // do zapisu danych
        do {
            lag_ss << b - a << ";" << std::endl;
            d_prev = d;
            a_sol.x = a;
            b_sol.x = b;
            c_sol.x = c;
            a_sol.fit_fun(ff, ud1);
            b_sol.fit_fun(ff, ud1);
            c_sol.fit_fun(ff, ud1);

            double l = m2d(a_sol.y) * (pow(b, 2) - pow(c, 2));
            l += m2d(b_sol.y) * (pow(c, 2) - pow(a, 2));
            l += m2d(c_sol.y) * (pow(a, 2) - pow(b, 2));
            double m = m2d(a_sol.y) * (b - c);
            m += m2d(b_sol.y) * (c - a);
            m += m2d(c_sol.y) * (a - b);

            if (m <= 0) {
                Xopt.flag = 0;
                break;
            }
            d = 0.5 * (l / m);
            d_sol.x = d;
            d_sol.fit_fun(ff, ud1);
            if (a < d && d < c) {
                if (d_sol.y < c_sol.y) {
                    //a = a;
                    b = c;
                    c = d;
                } else {
                    a = d;
                    //c = c;
                    //b = b;
                }
            } else if (c < d && d < b) {
                if (d_sol.y < c_sol.y) {
                    a = c;
                    c = d;
                    //b = b;
                } else {
                    //a = a;
                    //c = c;
                    b = d;
                }
            } else {
                Xopt.flag = 0;
                break;
            }
            i += 1;
            if (solution::f_calls > Nmax) {
                Xopt.flag = 0;
                break;
            }
        } while (!(b - a < epsilon or abs(d - d_prev) <= gamma));
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
    } catch (string ex_info) {
        throw ("solution lag(...):\n" + ex_info);
    }
}

solution HJ(matrix (*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax,
            matrix ud1, matrix ud2) {
    //etap probny -> badamy lokalne zachowanie funkcji wokol[ [unktu bazowego i uzywamy tego punktu bazowaego
    //etap roboczy -> przenosimys sie do nowego punktu i w nim robimy znowu etap probny

    try {
        solution Xopt;
        //Tu wpisz kod funkcji
        solution XB, X;
        XB.x = x0;
        XB.fit_fun(ff, ud1, ud2);
        std::stringstream HJ_ss; // do zapisu danych
        while (s > epsilon) {
            HJ_ss << m2d(XB.x(0)) << ";" <<  m2d(XB.x(1))<< ";\n";
            X = HJ_trial(ff, XB, s, ud1, ud2);
            if (X.y < XB.y) {
                do {
                    solution tempXB = XB;
                    XB = X;
                    X.x = 2.0 * XB.x - tempXB.x;
                    X.fit_fun(ff, ud1, ud2);
                    X = HJ_trial(ff, X, s, ud1, ud2);
                    if (solution::f_calls > Nmax) {
                        break;
                    }
                } while (X.y < XB.y);
                X = XB;
            } else {
                s = alpha * s;
            }
            if (solution::f_calls > Nmax) {
                Xopt.flag = 0;
                break;
            }
        }
        // zapis wynikow do pliku
        std::ofstream file(R"(C:\Users\Ania\CLionProjects\Optymalizacja\Optymalizacja\lab2-analiza\lab2-HJ-przyklad.txt)");
        if (file.is_open()) {
        	file << HJ_ss.str();
        	file.close();
        }else {
        	cerr << "Nie udało się otworzyć pliku do zapisu.\n";
        }
        return X;
    } catch (string ex_info) {
        throw ("solution HJ(...):\n" + ex_info);
    }

}

solution HJ_trial(matrix (*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2) {
    try {
        int n = get_len(XB.x);
        matrix e = matrix(n, n);
        for (int i = 0; i < n; i++) {
            e(i, i) = 1.0;
        }
        //Tu wpisz kod funkcji
        for (int i = 0; i < n; i++) {
            solution f1_sol, f2_sol;
            f1_sol.x = XB.x + s * e[i];
            f1_sol.fit_fun(ff, ud1, ud2);

            f2_sol.x = XB.x - s * e[i];
            f2_sol.fit_fun(ff, ud1, ud2);
            if (f1_sol.y < XB.y) {
                XB = f1_sol;
            } else if (f2_sol.y < XB.y) {
                XB = f2_sol;
            }

        }

        return XB;
    } catch (string ex_info) {
        throw ("solution HJ_trial(...):\n" + ex_info);
    }
}
double max(const matrix& m) {
    int len = get_len(m);  // Zakładam, że get_len() zwraca liczbę elementów w macierzy
    double result = 0.0;

    for (int i = 0; i < len; ++i) {
        result = std::max(result, std::abs(m(i)));
    }

    return result;
}
solution Rosen(matrix (*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon,
               int Nmax, matrix ud1, matrix ud2) {
    //wektor p to wektor porazek jak moje zycie
    //i nr czegos i J nr bazu? nwm do ktorego alg to bo nie sluchalam
    try {


        solution Xopt;
        //Tu wpisz kod funkcji
        int n = get_len(x0);
        std::stringstream Rosen_ss; // do zapisu danych
        matrix e = matrix(n, n);
        for (int i = 0; i < n; i++) {
            e(i, i) = 1.0;
        }

        int i = 0;
        matrix lambda = matrix(n, new double[n]{0.0});
        matrix p = matrix(n, new double[n]{0.0});
        matrix s = s0;

        solution XB;
        XB.x = x0;
        XB.fit_fun(ff, ud1, ud2);
        while (max(s) > epsilon) {
            Rosen_ss << m2d(XB.x(0)) << ";" <<  m2d(XB.x(1))<< ";\n";
            for (int j = 0; j < n; j++) {
                solution temp;
                temp.x = XB.x + s(j) * e[j];
                temp.fit_fun(ff, ud1, ud2);
                if (temp.y < XB.y) {
                    XB = temp;
                    lambda(j) += s(j);
                    s(j) *= alpha;
                } else {
                    s(j) *= -beta;
                    p(j) += 1;
                }
            }
            i++;
            Xopt = XB;
            bool changeDirection = true;
            for (int j = 0; j < n; j++) {
                if (lambda(j) == 0 || p(j) == 0) {
                    changeDirection = false;
                    break;
                }
            }
            if (changeDirection) {
                int l = 0;
                matrix lambdaM(n, n);
                for (int j = 0; j < n; j++) {
                    for (int k = 0; k <= j; k++) {
                        lambdaM(j, k) = lambda(l);
                    }
                    l++;
                }
                matrix Q = e * lambdaM;
                matrix V = matrix(n, 1);
                V = Q[0];
                e[0] = V/ norm(V);
                for (int j = 1; j < n; j++) {
                    matrix sum = matrix(n, new double[n]{0.0});
                    for (int k = 0; k<j; k++) {
                        sum = sum + (trans(Q[j]) * e[k]) * e[k];
                    }
                    V = Q[j] - sum;
                    e[j] = V/ norm(V);
                }
                lambda = matrix(n, new double[n]{0.0});
                p = matrix(n, new double[n]{0.0});
                s = s0;
            }
        }


        //macierz D to macierz kierunkow  d1 d2 i potem macierz
        //zlozona z lambd
        //zaimplementowac w 2D

        // zapis wynikow do pliku
        std::ofstream file(R"(C:\Users\Ania\CLionProjects\Optymalizacja\Optymalizacja\lab2-analiza\lab2-Rosen-przyklad.txt)"); //musialam dac cala sciezke bo nie dzialalo xd
        if (file.is_open()) {
            file << Rosen_ss.str();
            file.close();
        }else {
            cerr << "Nie udało się otworzyć pliku do zapisu.\n";
        }

        return Xopt;
    } catch (string ex_info) {
        throw ("solution Rosen(...):\n" + ex_info);
    }
}

solution pen(matrix (*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1,
             matrix ud2) {
    try {
        solution XB;
        XB.x = x0;
        XB.fit_fun(ff, ud1, c);

        solution XT;
        XT = XB;
        double s = ud2(0);
        double alpha = ud2(1);
        double beta = ud2(2);
        double gamma = ud2(3);
        double delta = ud2(4);

        while(true) {
            XT = sym_NM(ff, XB.x, s, alpha, beta, gamma, delta, epsilon, Nmax, ud1, c);
            c*=dc;
            if(solution::f_calls > Nmax) {
                XT.flag=0;
                throw std::string("Za duzo wykonan funkcji");
            }
            if (norm(XT.x - XB.x)<epsilon) break;
            XB=XT;
        }
        return XT;
    } catch (string ex_info) {
        throw ("solution pen(...):\n" + ex_info);
    }
}

double maxW(vector<solution>& sim, int i_min) {
    double result = 0.0;
    for (int i = 0; i<sim.size(); i++) {
        double normal = norm(sim[i_min].x - sim[i].x);
        result = max(result, normal);
    }
    return result;
}

solution sym_NM(matrix (*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma,
                double delta, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        int n = get_len(x0);
        matrix d = matrix(n, n);
        for (int i = 0; i < n; i++) {
            d(i, i) = 1.0;
        }


        vector<solution> sim;
        sim.resize(n+1);
        sim[0].x = x0;
        sim[0].fit_fun(ff, ud1, ud2);

        for (int i = 1; i < sim.size(); i++) {
            sim[i].x = sim[0].x + s*d[i-1];
            sim[i].fit_fun(ff, ud1, ud2);
        }
        int i_min{};
        int i_max{};
        while(maxW(sim, i_min)>=epsilon) {
            i_min =0;
            i_max =0;
            for (int j = 1; j < sim.size(); j++) {
                if(sim[j].y < sim[i_min].y) {
                    i_min = j;
                }
                if(sim[j].y > sim[i_max].y) {
                    i_max = j;
                }
            }
            matrix sim_SrC{};
            for (int i=0; i<sim.size(); i++) {
                if (i==i_max)
                    continue;
                sim_SrC = sim_SrC + sim[i].x;
            }
            sim_SrC = sim_SrC / sim.size();

            solution sim_reflected{};
            sim_reflected.x = sim_SrC + alpha * (sim_SrC - sim[i_max].x);
            sim_reflected.fit_fun(ff, ud1, ud2);

            if(sim_reflected.y < sim[i_max].y) {
                solution sim_expansion{};
                sim_expansion.x = sim_SrC + gamma * (sim_reflected.x - sim_SrC);
                sim_expansion.fit_fun(ff, ud1, ud2);
                if(sim_expansion.y < sim_reflected.y) {
                    sim[i_max] = sim_expansion;
                }
                else {
                    sim[i_max] = sim_reflected;
                }
            } else {
                if (sim[i_min].y <= sim_reflected.y && sim_reflected.y < sim[i_max].y ) {
                    sim[i_max] = sim_reflected;
                } else {
                    solution sim_narrowed{};
                    sim_narrowed.x = sim_SrC + beta * (sim[i_max].x - sim_SrC);
                    sim_narrowed.fit_fun(ff, ud1, ud2);
                    if (sim_narrowed.y >= sim[i_max].y) {
                        for (int i = 0; i<sim.size(); i++) {
                            if (i==i_min)
                                continue;
                            sim[i].x = delta * (sim[i].x + sim[i_min].x);
                            sim[i].fit_fun(ff, ud1, ud2);
                        }
                    }else {
                        sim[i_max] = sim_narrowed;
                    }
                }
            }
            if(solution::f_calls > Nmax) {
                sim[i_min].flag=0;
                throw std::string("Za duzo wykonan funkcji");
            }

         }
        return sim[i_min];
    } catch (string ex_info) {
        throw ("solution sym_NM(...):\n" + ex_info);
    }
}

solution SD(matrix (*ff)(matrix, matrix, matrix), matrix (*gf)(matrix, matrix, matrix), matrix x0, double h0,
            double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        std::stringstream test_ss;
        solution Xopt, D;
        Xopt.x = x0;
        D = Xopt;
        matrix d(2, 1);
        while(true) {
            test_ss << Xopt.x(0) << ";" << Xopt.x(1) << ";\n";
            Xopt.grad(gf, ud1, ud2);
            d = -Xopt.g;
            if (h0<=0) {
                matrix h_data(2, 2);
                h_data.set_col(Xopt.x, 0);
                h_data.set_col(d, 1);
                solution H = golden(ff, 0, 1, epsilon, Nmax, ud1, h_data);
                solution::f_calls = 0;
                D.x = D.x + H.x * d;
            }else
                D.x = D.x + h0 * d;
            if(solution::g_calls > Nmax) {
                // std::ofstream file3(R"(C:\Users\aneta\CLionProjects\Optymalizacja\Optymalizacja\lab4-analiza\lab4-zkSD.txt)");
                // if (file3.is_open()) {
                //     file3 << test_ss.str();
                //     file3.close();
                // }else {
                //     cerr << "Nie udało się otworzyć pliku do zapisu.\n";
                // }
                D.fit_fun(ff, ud1, ud2);
                return D;
                throw std::string("Za duzo wykonan funkcji");
            }
            if (norm(D.x - Xopt.x) <= epsilon)
                break;
            Xopt = D;
        }
        // std::ofstream file3(R"(C:\Users\aneta\CLionProjects\Optymalizacja\Optymalizacja\lab4-analiza\lab4-zkSD.txt)");
        // std::ofstream file3("C:\\Users\\Animatt\\CLionProjects\\Optymalizacja\\Optymalizacja\\lab4-analiza\\lab4-05SD.txt"); //musialam dac cala sciezke bo nie dzialalo xd
        // std::ofstream file3("C:\\Users\\Animatt\\CLionProjects\\Optymalizacja\\Optymalizacja\\lab4-analiza\\lab4-12SD.txt"); //musialam dac cala sciezke bo nie dzialalo xd
        // std::ofstream file3("C:\\Users\\Animatt\\CLionProjects\\Optymalizacja\\Optymalizacja\\lab4-analiza\\lab4-zkSD.txt"); //musialam dac cala sciezke bo nie dzialalo xd
        // if (file3.is_open()) {
        //     file3 << test_ss.str();
        //     file3.close();
        // }else {
        //     cerr << "Nie udało się otworzyć pliku do zapisu.\n";
        // }
        D.fit_fun(ff, ud1, ud2);
        return D;
    } catch (string ex_info) {
        throw ("solution SD(...):\n" + ex_info);
    }
}

solution CG(matrix (*ff)(matrix, matrix, matrix), matrix (*gf)(matrix, matrix, matrix), matrix x0, double h0,
            double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        std::stringstream test_ss;
        solution Xopt, Xp;
        Xopt = x0;
        Xopt.grad(gf, ud1, ud2);
        matrix d = -Xopt.g;
        matrix di;
        do {
            test_ss << Xopt.x(0) << ";" << Xopt.x(1) << ";\n";
            Xp = Xopt;
            d = di;
            Xopt.grad(gf, ud1, ud2);
            double beta = pow(norm(Xopt.g),2) / pow(norm(Xp.g),2);
            di = -Xopt.g + beta * d;
            if (h0<=0) {
                matrix h_data(2, 2);
                h_data.set_col(Xp.x, 0);
                h_data.set_col(di, 1);
                solution H = golden(ff, 0, 1, epsilon, Nmax, ud1, h_data);
                solution::f_calls = 0;
                matrix h = H.x;
                Xopt.x = Xp.x + h * di;
            }else {
                Xopt.x = Xp.x + h0 * di;
            }
            if(solution::g_calls > Nmax) {
                Xopt.fit_fun(ff, ud1, ud2);
                return Xopt;
                throw std::string("Za duzo wykonan funkcji");
            }

        }
        while(norm(Xopt.x - Xp.x) > epsilon);
        std::ofstream file3(R"(C:\Users\aneta\CLionProjects\Optymalizacja\Optymalizacja\lab4-analiza\lab4-zkCG.txt)");
        // std::ofstream file3("C:\\Users\\Animatt\\CLionProjects\\Optymalizacja\\Optymalizacja\\lab4-analiza\\lab4-zkCG.txt"); //musialam dac cala sciezke bo nie dzialalo xd
        if (file3.is_open()) {
            file3 << test_ss.str();
            file3.close();
        }else {
            cerr << "Nie udało się otworzyć pliku do zapisu.\n";
        }
        Xopt.fit_fun(ff, ud1, ud2);
        return Xopt;
    } catch (string ex_info) {
        throw ("solution CG(...):\n" + ex_info);
    }
}

solution Newton(matrix (*ff)(matrix, matrix, matrix), matrix (*gf)(matrix, matrix, matrix),
                matrix (*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1,
                matrix ud2) {
    try {
        std::stringstream test_ss;
        solution Xopt, D;
        D.x = x0;
        Xopt = D;
        matrix d;
        while(true) {
            test_ss << Xopt.x(0) << ";" << Xopt.x(1) << ";\n";
            D.grad(gf, ud1, ud2);
            D.hess(Hf, ud1, ud2);
            d = -inv(D.H) * D.g;
            if (h0<=0) {
                matrix h_data(2, 2);
                h_data.set_col(D.x, 0);
                h_data.set_col(d, 1);
                solution H = golden(ff, 0, 1, epsilon, Nmax, ud1, h_data);
                solution::f_calls = 0;
                Xopt.x = D.x + H.x * d;
            }else
                Xopt.x = D.x + h0 * d;
            if(solution::g_calls > Nmax) {
                throw std::string("Za duzo wykonan funkcji");
            }
            if(norm(Xopt.x - D.x) <= epsilon) {
                break;
            }
            D = Xopt;

        }
        std::ofstream file3(R"(C:\Users\aneta\CLionProjects\Optymalizacja\Optymalizacja\lab4-analiza\lab4-zkNewton.txt)");
        // std::ofstream file3("C:\\Users\\Animatt\\CLionProjects\\Optymalizacja\\Optymalizacja\\lab4-analiza\\lab4-zkNewton.txt"); //musialam dac cala sciezke bo nie dzialalo xd
        if (file3.is_open()) {
            file3 << test_ss.str();
            file3.close();
        }else {
            cerr << "Nie udało się otworzyć pliku do zapisu.\n";
        }
        Xopt.fit_fun(ff, ud1, ud2);
        return Xopt;
    } catch (string ex_info) {
        throw ("solution Newton(...):\n" + ex_info);
    }
}

solution golden(matrix (*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1,
                matrix ud2) {
    try {
        solution Xopt;
        double alpha = (pow(5,0.5)-1) / 2;
        double a0 = a;
        double b0 = b;
        double c0 = b0 - alpha*(b0 - a0);
        double d0 = a0 + alpha*(b0 - a0);
        do {
            solution fc, fd;
            fc.x = c0;
            fd.x = d0;
            fc.fit_fun(ff, ud1, ud2);
            fd.fit_fun(ff, ud1, ud2);

            if(fc.y < fd.y) {
                b0 = d0;
                d0 = c0;
                c0 = b0 - alpha*(b0 - a0);
            }else {
                a0 = c0;
                c0 = d0;
                d0 = a0 + alpha*(b0 - a0);
            }
            if(solution::f_calls > Nmax) {
                throw std::string("Za duzo wykonan funkcji");
            }

        }
        while(b0-a0>=epsilon);

        Xopt.x = (a0 + b0)/2;
        return Xopt;
    } catch (string ex_info) {
        throw ("solution golden(...):\n" + ex_info);
    }
}

solution Powell(matrix (*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    } catch (string ex_info) {
        throw ("solution Powell(...):\n" + ex_info);
    }
}

solution EA(matrix (*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0,
            double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    } catch (string ex_info) {
        throw ("solution EA(...):\n" + ex_info);
    }
}
